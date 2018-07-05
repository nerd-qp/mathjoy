
#ifndef _FINITE_ELEMENT_METHOD_
#define _FINITE_ELEMENT_METHOD_
/*
  Include some dependent libs
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/range/irange.hpp>
#include <gsl/gsl_integration.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <iostream>
#include <type_traits>
#include <algorithm>
#include <numeric>

#include "utils/utils.hpp"

/*
  Some macros
 */
#define DEFAULT_INTERPOLATION_NODES_NUM 5
#define COMPUTE_ALGEBRAIC_PRECISION(X) (2 * (X) - 1)
#define NORMALIZED_RANGE_LEFT 0
#define NORMALIZED_RANGE_RIGHT 1
#define DEFAULT_DEGREE 1


/*
  This class is
  to solve the differential equation of the form -d(p du) + q du = f
  With finite element method(linear)
*/

namespace linear {
  // data structure to hold a(phi_i, phi_j)
  class finite_element_solver {
  private:
    // private fields

    // _stiffness_matrix is a sparse matrix with only half of the nodes
    // represented by a vector, with operator(i, j) to act like a matrix
    std::vector<double> _stiffness_matrix;    // given the symmetric, target
    std::vector<double> _right_f_phi;         // target
    std::vector<double> _interpolation_nodes; // has one more element than hi
    gsl_integration_glfixed_table *_integral_table; // The gsl table

    // Most of these private fields are set in the constructor
    std::vector<double> _unit_interval;             // inteval of each unit

    // Three functions
    double (*_p)(double);
    double (*_q)(double);
    double (*_f)(double);

    // The precision
    int _precision_gq; // The precision of Gaussian quadrature
    int _n;            // the size of hi
    int _a, _b;        // ranges [a, b]
    bool _state_stiff_mat;   // Default to false at first
    bool _state_right_f_phi; // Default to false at first

    // private methods, this method assumes the stiffness matrix is symmetric
    // the core of the class, compute the stiffness matrix and the (f, phi)
    bool compute_stiffness_matrix();
    bool compute_right_f_phi();

  public:
    finite_element_solver();

    // the basic input needed to construct
    // stiffness matrix and right side (f, phi)
    finite_element_solver(
                          // The range of the equation
                          const double &a, const double &b,
                          // the unit intevals
                          const std::vector<double> &hi,
                          // The three functions
                          double (*p)(double),
                          double (*q)(double),
                          double (*f)(double),
                          // precision for Gaussian quadrature
                          int n = DEFAULT_INTERPOLATION_NODES_NUM
                          );

    // The structure to use with gsl integral
    struct _p_func_params {
      double (*_p_func)(double);
      double x;
      double h;
    };

    // carry one more function with it
    struct _q_func_params : _p_func_params {
      double (*_q_func)(double);
    };

    /*
      The integral function that takes an x and h
      to map a function ranged in [x, x+h] -> [0, 1]
      return value of the compound function with parameter xi
    */
    static double p_func_integral(double xi, void *params);
    static double q_func_integral(double xi, void *params);

    // takes a vector of hi, which is the unit slice's length of i
    // returns the xi as vector
    // does check on boundaries
    static bool x_i_array(double a, double b,
                          const std::vector<double> &hi,
                          std::vector<double> &res);

    // accessor range from [0~n-1, 0~n-1]
    double& operator() (int, int);

    boost::numeric::ublas::matrix<double> get_stiffness_matrix();
    boost::numeric::ublas::vector<double> get_right_f_phi();
    boost::numeric::ublas::vector<double> get_solution();

  };
}



namespace general {
  // class FEM_solver {
  // private:
  //   // Computes the local stiffness matrix
  //   // No exception when index exceeds boundary
  //   virtual std::vector<double> _local_stiffness_compute(int index) = 0;

  //   // Right hand side
  //   // std::vector<double> right_f_phi();
  //   virtual void _right_f_phi() = 0;

  //   // Assemble total stiffness matrix
  //   virtual void _compute_stiffness_matrix() = 0;

  //   virtual void _solve() = 0;

  // public:
  //   virtual boost::numeric::ublas::vector<double> get_solution() = 0;

  // };

  // class FEM_solver_impl : public FEM_solver {
  // private:
  //   int _n;

  //   // The degree of the intepolation
  //   int _degree;

  // public:
  //   FEM_solver_impl(int n, int degree) : _n(n), _degree(degree) {}
  // };

  class FEM_solver {
  private:
    int _freedom;
    int _degree;
    double _left;
    double _right;

    std::function<double (double)> _p;
    std::function<double (double)> _q;
    std::function<double (double)> _f;
    // double (*_p)(double);
    // double (*_q)(double);
    // double (*_f)(double);

    // Most of these private fields are set in the constructor
    std::vector<double> _unit_interval;             // inteval of each unit
    std::vector<double> _interpolation_nodes;

    std::vector<double> _local_node;
    // Denominator of index k
    std::vector<double> _denominator;
    std::vector<std::function<double (double)>> _N;
    std::vector<std::function<double (double)>> _N_diff;

    int _precision_gq; // The precision of Gaussian quadrature
    int _n;            // the size of hi
    gsl_integration_glfixed_table *_integral_table; // The gsl table

    // TODO: speed would increase if we know before hand the nonzero size
    boost::numeric::ublas::mapped_matrix<double> _stiffness_matrix;
    boost::numeric::ublas::vector<double> _right_f_phi;

    // given index of the small region -> index
    // and index of local intepolation -> x
    inline int _trans_local_x(int index, int x);
    inline int _trans_local_y(int index, int y);

    inline void _compute_denominator_all();
    inline double _compute_denominator_at(int index);
    inline void _compute_canonical_nodes();

    // Compute and storage all LagrangePolynomial and their differentiation
    inline void _compute_LagrangePolynomial_all();

    // TODO: Need to take care to avoid intepolation point as negative
    // Or degree as negative, or 0

    // compute local stiffness given the starting and ending degree
    inline void _compute_local_stiffness(int index, int d_start, int d_end);

    // simple wrapper with default argument, since using field variable
    // that may not be initialized
    inline void _compute_local_stiffness(int index);

    inline void _compute_stiffness_matrix();

    inline void _compute_local_right_f_phi(int index, int d_start, int d_end);

    inline void _compute_local_right_f_phi(int index);

    inline void _compute_right_f_phi();


    /*
      Generate Lagrange Polynomial

      precondition:
      1. canonical nodes must be computed first
      2. the range, degree of the problem must be initialized

      purpose:
      1. provide index and return function

      warning: No checking for the function
    */
    std::function<double (double)>
    LagrangePolynomial(int index);

    // short for differentiation of LagrangePolynomial
    std::function<double (double)>
    LagrangePolynomial_diff(int index);

    // function parameters for integral
    struct func_params {
      std::function<double (double)> func;
      std::function<double (double)> prod_func;
      double x;
      double h;
    };

    // simple function for integral
    static double
    func_integral(double xi, void *params);

  public:
    //static double
    FEM_solver();
    // The default constructor should not be used
    FEM_solver(double left, double right,
               const std::vector<double> &hi,
               std::function<double (double)> p,
               std::function<double (double)> q,
               std::function<double (double)> f,
               int degree = DEFAULT_DEGREE,
               int precision = DEFAULT_INTERPOLATION_NODES_NUM);

    static bool
    x_i_array(double a, double b,
              const std::vector<double> &hi,
              std::vector<double> &res);


    boost::numeric::ublas::mapped_matrix<double>
    get_stiffness_matrix();

    boost::numeric::ublas::vector<double>
    get_right_f_phi();

    boost::numeric::ublas::vector<double>
    get_solution();

  };

}

#endif
