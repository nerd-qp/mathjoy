#include "finite_element_method.hpp"
#include <boost/numeric/ublas/io.hpp>
#include "matrix_solver/cholesky.hpp"
#include "solver/iterative_solve.hpp"
#include <cmath>

namespace general
{




  int FEM_solver::_trans_local_x(int index, int x) {
    return _degree * index + x - 1;
  }

  int FEM_solver::_trans_local_y(int index, int y) {
    return _degree * index + y - 1;
  }

  void FEM_solver::_compute_LagrangePolynomial_all() {
    for ( auto i : boost::irange(0, _degree + 1) ) {
      _N[i] = LagrangePolynomial(i);
      _N_diff[i] = LagrangePolynomial_diff(i);
    }
  }

  void FEM_solver::_compute_denominator_all() {
    for ( const auto & i : boost::irange(0, (int)_denominator.size()) )
      _denominator[i] = _compute_denominator_at(i);
  }

  double FEM_solver::_compute_denominator_at(int index) {
    double denominator = 1;

    auto start = _local_node.cbegin();
    // provided that the index doesn't exceed
    // TODO: can add exception here
    auto mid = _local_node.cbegin() + index;
    auto end = _local_node.cend();
    auto index_val = *mid;

    std::for_each(start, mid,
                  [&] ( double x ) {
                    denominator *= (index_val - x);
                  });

    std::for_each(mid+1, end,
                  [&] ( double x ) {
                    denominator *= (index_val - x);
                  });

    return denominator;
  }

  void FEM_solver::_compute_canonical_nodes() {
    // If canonical range is [0, 1]
    // Compute local nodes as 0, 1/p, 2/p, ..., 1
    //_local_node = std::vector<double>(_degree + 1);
    double steps = ((NORMALIZED_RANGE_RIGHT - NORMALIZED_RANGE_LEFT) /
                    static_cast<double>(_degree));

    for ( size_t i = 0; i < _local_node.size(); ++ i )
      _local_node[i] = NORMALIZED_RANGE_LEFT + steps * i;
  }

  void FEM_solver::_compute_local_stiffness(int index) {
    _compute_local_stiffness(index, 0, _degree+1);
  }

  void FEM_solver::_compute_local_stiffness(int index,
                                            int d_start,
                                            int d_end) {

    func_params params;
    params.func = _p;
    //auto LagrangePolynomial_diff(index);
    params.x = _interpolation_nodes[index];
    params.h = _unit_interval[index];

    gsl_function F;
    F.function = func_integral;

    for ( int i = d_start; i < d_end; ++ i ) {

      for ( int j = d_start; j < d_end; ++ j ) {

        const auto N_diff_i = _N_diff.cbegin() + i;
        const auto N_diff_j = _N_diff.cbegin() + j;

        // have to use const_iterator in order to captured by reference
        params.prod_func = [&](double xi) {
                             return (*N_diff_i)(xi) * (*N_diff_j)(xi);
                           };

        F.params = &params;

        _stiffness_matrix
          (_trans_local_x(index, i),
           _trans_local_y(index, j)) +=
          gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                                  NORMALIZED_RANGE_RIGHT,
                                  _integral_table) / params.h;

      }

    }

    params.func = _q;

    for ( int i = d_start; i < d_end; ++ i ) {

      for ( int j = d_start; j < d_end; ++ j ) {

        const auto N_i = _N.cbegin() + i;
        const auto N_j = _N.cbegin() + j;

        params.prod_func = [&](double xi) {
                             return (*N_i)(xi) * (*N_j)(xi);
                           };

        F.params = &params;

        _stiffness_matrix
          (_trans_local_x(index, i),
           _trans_local_y(index, j)) +=
          gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                                  NORMALIZED_RANGE_RIGHT,
                                  _integral_table) * params.h;
      }
    }

  }

  void FEM_solver::_compute_stiffness_matrix() {

    // Treating the first LPolynomial of index 0 special
    // Given precondition of u0
    int index = 0;
    _compute_local_stiffness(0, 1, _degree+1);

    for ( int i = 1; i < _unit_interval.size(); ++i ) {

      _compute_local_stiffness(i);

    }

  }

  void FEM_solver::_compute_local_right_f_phi(int index) {
    _compute_local_right_f_phi(index, 0, _degree+1);
  }

  void FEM_solver::_compute_local_right_f_phi(int index, int d_start, int d_end) {

    func_params params;
    params.func = _f;
    //auto LagrangePolynomial_diff(index);
    params.x = _interpolation_nodes[index];
    params.h = _unit_interval[index];

    gsl_function F;
    F.function = func_integral;

    for ( const auto &i : boost::irange(d_start, d_end) ) {

      const auto N_i = _N.cbegin() + i;

      params.prod_func = [&](double xi) {
                           return (*N_i)(xi);
                         };

      F.params = &params;

      _right_f_phi(_trans_local_x(index, i)) +=
        gsl_integration_glfixed(&F, NORMALIZED_RANGE_LEFT,
                                NORMALIZED_RANGE_RIGHT,
                                _integral_table) * params.h;

    }

  }

  void FEM_solver::_compute_right_f_phi() {

    int index = 0;
    _compute_local_right_f_phi(0, 1, _degree+1);

    for ( const auto &i : boost::irange(1, (int)_unit_interval.size() ) ) {

      _compute_local_right_f_phi(i);

    }
  }

  std::function<double (double)> FEM_solver::LagrangePolynomial(int index) {
    auto start = _local_node.cbegin();
    // provided that the index doesn't exceed
    // TODO: can add exception here
    auto mid = _local_node.cbegin() + index;
    auto end = _local_node.cend();

    // Do we really have closure?
    // Very impressive work to include lambda in cpp
    return ([=] ( double x ) {
              double res = 1;
              std::for_each(start, mid, [&] ( double y ) {
                                          res *= (x - y);
                                        });

              std::for_each(mid+1, end, [&] ( double y ) {
                                          res *= (x - y);
                                        });
              return res / _denominator[index];
            });
  }

  std::function<double (double)> FEM_solver::LagrangePolynomial_diff(int index) {
    auto start = _local_node.cbegin();
    auto mid = _local_node.cbegin() + index;
    auto end = _local_node.cend();

    // First nasty version
    return ([=] ( double x ) {
              double res = 0;
              for ( auto i = start; i != end; ++ i ) {
                if ( i != mid ) {
                  double part = 1;
                  for ( auto j = start; j != end; ++ j ) {
                    if (j != i && j != mid ) {
                      part *= ( x - *j );
                    }
                  }
                  res += part;
                }
              }
              return res / _denominator[index];
            });

    // second better version
  }

  double FEM_solver::func_integral(double xi, void *params) {
    func_params *the_fun = static_cast<func_params*>(params);
    return (the_fun->func(xi * the_fun->h + the_fun->x) * the_fun->prod_func(xi));
  }

  FEM_solver::FEM_solver() {};

  FEM_solver::FEM_solver(double left, double right,
                         const std::vector<double> &hi,
                         std::function<double (double)> p,
                         std::function<double (double)> q,
                         std::function<double (double)> f,
                         int degree,
                         int precision)
    : _left(left), _right(right),

      // Initialize with size to avoid performance issue with push_back
      _degree(degree), _local_node(degree+1),

      // default to zero to avoid initialization by hand
      _stiffness_matrix(hi.size() * degree, hi.size() * degree),

      _right_f_phi(hi.size() * degree, 0),

      // Initialize with size to avoid performance issue with push_back
      _denominator(degree+1), _p(p), _q(q), _f(f), _n(hi.size()),

      // precision are to be computed to suit GQ integral
      _precision_gq(COMPUTE_ALGEBRAIC_PRECISION(precision)),

      // integral table pointer defaults to a null pointer
      _integral_table(NULL), _N(degree+1), _N_diff(degree+1),

      // Save inteval information
      _unit_interval(hi)
  {

    _compute_canonical_nodes();
    _compute_denominator_all();
    _compute_LagrangePolynomial_all();
    x_i_array(left, right, hi, _interpolation_nodes);

    // allocate for integral computati
    _integral_table =
      gsl_integration_glfixed_table_alloc(_precision_gq);

    _compute_stiffness_matrix();
    _compute_right_f_phi();

  };

  bool FEM_solver::x_i_array(double a, double b,
                             const std::vector<double> &hi,
                             std::vector<double> &res) {

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

  using namespace boost::numeric::ublas;

  mapped_matrix<double>
  FEM_solver::get_stiffness_matrix() {
    return _stiffness_matrix;
  }

  vector<double>
  FEM_solver::get_right_f_phi() {
    return _right_f_phi;
  }

  vector<double>
  FEM_solver::get_solution() {
    auto my_matrix = this->get_stiffness_matrix();
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
