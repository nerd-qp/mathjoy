#include "finite_element_method.hpp"
#include <boost/numeric/ublas/io.hpp>
#include "matrix_solver/cholesky.hpp"
#include "solver/iterative_solve.hpp"
#include <cmath>

namespace linear
{

  // This function assumes the variables are already initialized
  bool finite_element_solver::compute_stiffness_matrix() {
    // Initialize func_params with p and q
    _p_func_params p;
    _q_func_params q;
    p._p_func = _p;
    q._p_func = _q; // q carries one more function, depending on (i,j)
    double h, x;

    // edge case for (n, n) => half of the area

    // set h, x
    h = _unit_interval[_n-1];       // h(n)
    x = _interpolation_nodes[_n-1]; // x(n-1)
    p.h = q.h = h;
    p.x = q.x = x;

    // set q additional function
    q._q_func = +[](double x) { return x*x; };

    // Set F for p
    gsl_function F;
    F.function = p_func_integral;
    F.params = &p;
    (*this)(_n-1, _n-1) +=
      gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                              NORMALIZED_RANGE_RIGHT, _integral_table) / p.h;

    // Set F for q
    F.function = q_func_integral;
    F.params = &q;
    (*this)(_n-1, _n-1) +=
      gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                              NORMALIZED_RANGE_RIGHT, _integral_table) * q.h;

    for ( int i : boost::irange(0, _n-1) ) {
      // for (i, i) x(j-1) to x(j)
      q.h = p.h = _unit_interval[i];
      q.x = p.x = _interpolation_nodes[i];
      F.function = p_func_integral;
      F.params = &p;
      (*this)(i, i) +=
        gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                                NORMALIZED_RANGE_RIGHT, _integral_table) / p.h;

      q._q_func = +[](double x) { return x*x; };
      F.function = q_func_integral;
      F.params = &q;
      (*this)(i, i) +=
        gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                                NORMALIZED_RANGE_RIGHT, _integral_table) * q.h;

      q.h = p.h = _unit_interval[i+1];
      q.x = p.x = _interpolation_nodes[i+1];
      F.function = p_func_integral;
      F.params = &p;
      (*this)(i, i) +=
        gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                                NORMALIZED_RANGE_RIGHT, _integral_table) / p.h;

      q._q_func = +[](double x) { return (1-x)*(1-x); };
      F.function = q_func_integral;
      F.params = &q;
      (*this)(i, i) +=
        gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                                NORMALIZED_RANGE_RIGHT, _integral_table) * q.h;

      // for (i, i + 1)
      // q.h and p.h, q.x and p.x are the same as above
      q._q_func = +[](double x) { return (1-x)*x; };
      F.params = &q;
      (*this)(i, i+1) +=
        gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                                NORMALIZED_RANGE_RIGHT, _integral_table) * q.h;

      F.function = p_func_integral;
      F.params = &p;
      (*this)(i, i+1) -=
        gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                                NORMALIZED_RANGE_RIGHT, _integral_table) / p.h;
    }

    return true;
  }

  bool finite_element_solver::compute_right_f_phi() {
     // _f_func_params would be similar to q, so using _q_func_params
    _q_func_params f;
    gsl_function F;
    F.function = q_func_integral;

    // setting edge case for (n,n)
    f.h = _unit_interval[_n-1];
    f.x = _interpolation_nodes[_n-1];
    f._p_func = _f;
    f._q_func = +[](double x) { return x; };
    F.params = &f;

    _right_f_phi[_n-1] +=
      gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                              NORMALIZED_RANGE_RIGHT, _integral_table) * f.h;

    // loop for 1 -> n-1
    for ( auto i : boost::irange(0, _n-1) ) {

      // integral from x(i-1) -> x(i)
      f.h = _unit_interval[i];
      f.x = _interpolation_nodes[i];
      f._q_func = +[](double x) { return x; };
      F.params = &f;

      _right_f_phi[i] +=
        gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                                NORMALIZED_RANGE_RIGHT, _integral_table) * f.h;

      // integral from x(i) -> x(i+1)
      f.h = _unit_interval[i+1];
      f.x = _interpolation_nodes[i+1];
      f._q_func = +[](double x) { return (1-x); };
      F.params = &f;

      _right_f_phi[i] +=
        gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                                NORMALIZED_RANGE_RIGHT, _integral_table) * f.h;
    }

    return true;
  }

  bool finite_element_solver::x_i_array(double a, double b,
                                        const std::vector<double> &hi,
                                        std::vector<double> &res)
  {

    double tmp = a;
    for (const auto &h : hi) {
      if ( tmp <= b )
        res.emplace_back(tmp);
      else
        return false;
      tmp += h;
    }
    res.emplace_back(tmp);

    return utils::almost_equal(tmp, b);
  }

  double finite_element_solver::p_func_integral(double xi, void * params) {
    _p_func_params * the_fun = static_cast<_p_func_params*>(params);
    return the_fun->_p_func(xi*the_fun->h + the_fun->x);
  }

  double finite_element_solver::q_func_integral(double xi, void * params) {
    _q_func_params * the_fun = static_cast<_q_func_params*>(params);
    return the_fun->_p_func(xi*the_fun->h + the_fun->x) * the_fun->_q_func(xi);
  }

  finite_element_solver::finite_element_solver()
    : _precision_gq(0), _integral_table(NULL),
      _p(NULL), _q(NULL), _f(NULL), _n(0),
      _state_stiff_mat(false), _state_right_f_phi(false) {}

  finite_element_solver::finite_element_solver
  (const double &a, const double &b,
   const std::vector<double> &hi,
   double (*p)(double),
   double (*q)(double),
   double (*f)(double),
   int n) : _a(a), _b(b), _p(p), _q(q), _f(f),
            _precision_gq(COMPUTE_ALGEBRAIC_PRECISION(n)),
            _unit_interval(hi), _n(hi.size()),
            _stiffness_matrix(std::vector<double>((hi.size()+1)
                                                  *hi.size()/2, 0)),
            _right_f_phi(std::vector<double>(hi.size(), 0)),
            _state_stiff_mat(false),
            _state_right_f_phi(false){}

  double& finite_element_solver::operator() (int i, int j) {
    if ( i > j ) return (*this)(j, i);
    // mapping from (i,j) -> symtric matrix
    else return _stiffness_matrix[((2*_n - i + 1) * i) / 2 + j - i];
  }

  boost::numeric::ublas::matrix<double>
  finite_element_solver::get_stiffness_matrix() {
    using namespace boost;
    using namespace numeric::ublas;
    /*
      Calculate the stiffness matrix if not already
     */
    if ( !_state_stiff_mat ) {
      if ( !_state_right_f_phi ) {
        // avoid reinitialization, need a better way to do this though
        _integral_table =
          gsl_integration_glfixed_table_alloc(_precision_gq);

        x_i_array(_a, _b, _unit_interval, _interpolation_nodes);
      }
      compute_stiffness_matrix();
      _state_stiff_mat = true;
    }

    /*
      shuffle the stiffness matrix and return
     */
    matrix<double> stiffness_matrix (_n, _n);
    for ( auto i : irange(0, _n) )
      for ( auto j : irange(0, _n))
        stiffness_matrix (i, j) = (*this)(i, j);

    // compute
    return stiffness_matrix;
  }

  boost::numeric::ublas::vector<double>
  finite_element_solver::get_right_f_phi() {
    using namespace boost;
    using namespace numeric::ublas;
    if ( !_state_right_f_phi ) {
      if ( !_state_stiff_mat ) {
        // avoid reinitialization, need a better way to do this though
        _integral_table =
          gsl_integration_glfixed_table_alloc(_precision_gq);
        x_i_array(_a, _b, _unit_interval, _interpolation_nodes);
      }
      compute_right_f_phi();
      _state_right_f_phi = true;
    }

    vector<double> right_f_phi(_right_f_phi.size());
    for ( auto i : irange(0, (int)right_f_phi.size()) )
      right_f_phi(i) = _right_f_phi[i];
    return right_f_phi;
  }

  boost::numeric::ublas::vector<double>
  finite_element_solver::get_solution() {
    auto my_matrix =  this->get_stiffness_matrix();
    auto my_vector = this->get_right_f_phi();

    boost::numeric::ublas::permutation_matrix<size_t> pm(my_matrix.size1());
    lu_factorize(my_matrix, pm);
    lu_substitute(my_matrix, pm, my_vector);

    // incomplete_cholesky_decompose(my_matrix);
    // cholesky_solve(my_matrix, my_vector, ublas::lower());
    return my_vector;

    // auto x = my_vector;
    // []Gauss_Seidel_Method(my_matrix, my_vector, x);

    // auto x = my_vector;
    // Conjugate_Gradient_Method(my_matrix, my_vector, x);

    //return x;
  }

}
