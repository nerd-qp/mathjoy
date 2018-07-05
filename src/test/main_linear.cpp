

//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
// gsl libraries
//#include <gsl/gsl_integration.h>

/*
  Include libraries
*/

// std libraries
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>

// boost libraries
#include <boost/range/irange.hpp>

// self
#include "finite_element_method.hpp"
#include "testing_functions.hpp"

// for interpolation nodes
#define MAXIMUM_NODES 100

// int main() {
//   using namespace general;
//   FEM_solver x(1,1);
// }

int main() {
  int i = 10; double b = 1; double a = 0;
  using namespace linear;

  std::vector<double> hi(i, (b-a)/(double)i), xi;

  auto pre = std::chrono::high_resolution_clock::now();
  finite_element_solver solver(a, b, hi, &test_p, &test_q, &test_f);
  std::cout << "Time used to initialize: " <<
    std::chrono::duration_cast<std::chrono::duration<double>>
    (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;

  std::cout << hi[0] * solver.get_stiffness_matrix() << std::endl;
  std::cout << "Total time used: " <<
    std::chrono::duration_cast<std::chrono::duration<double>>
    (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;
  // int a = 0, b = 1;
  // //for ( int i = 1; i <= MAXIMUM_NODES; i*=10 ) {
  // for ( int i = 1; i <= MAXIMUM_NODES; i = (i*2 + 1) ) {
  //   //std::vector<double> hi(100, 1e-2), xi;
  //   //std::vector<double> hi(10, 1e-1), xi;
  //   std::vector<double> hi(i, (b-a)/(double)i), xi;

  //   using namespace linear;
  //   finite_element_solver::x_i_array(a, b, hi, xi);

  //   std::cout << "Interpolation nodes num: " << i << std::endl;

  //   // compute and time it
  //   auto pre = std::chrono::high_resolution_clock::now();
  //   finite_element_solver solver(a, b, hi, test_p, test_q, test_f, 100);
  //   boost::numeric::ublas::vector<double> my_vector = solver.get_solution();
  //   std::cout << "Time used: " <<
  //     std::chrono::duration_cast<std::chrono::duration<double>>
  //     (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;

  //   double total_error = 0;
  //   //std::cout << my_vector << std::endl;

  //   for ( auto i : boost::irange(1, (int)xi.size()) )
  //     total_error += fabs(result_f(xi[i]) - my_vector(i-1));

  //   std::cout << "Total error is: " << total_error << std::endl;
  //   std::cout << "Mean error is: "
  //     // << total_error / hi.size()
  //             << total_error / hi.size() << "\n"
  //             << total_error / (sqrt(hi[0]))
  //             << "\n" << std::endl;
  // }
}

// testing area
//finite_element_solver::x_i_array(a, b, hi, xi);
//boost::numeric::ublas::matrix<double> my_matrix =  solver.get_stiffness_matrix();
//boost::numeric::ublas::vector<double> my_vector = solver.get_right_f_phi();

//boost::numeric::ublas::vector<double> sol;
//boost::numeric::ublas::permutation_matrix<size_t> pm(my_matrix.size1());
//lu_factorize(my_matrix, pm);
//lu_substitute(my_matrix, pm, my_vector);

// incomplete_cholesky_decompose(my_matrix);
// cholesky_solve(my_matrix, my_vector, ublas::lower());
  // func_params tmp_fun;
  // tmp_fun.func = test_f;
  // tmp_fun.x = 3;
  // tmp_fun.h = 10;
  // gsl_function F;
  // F.function = &func_integral;
  // F.params = &tmp_fun;



  // // Do note that this is the most time comsumming part
  // // Please check out Gauss-Legendre about how to compute integral
  // int n = 10000;
  // auto pre = std::chrono::high_resolution_clock::now();
  // gsl_integration_glfixed_table *w_table =
  //   gsl_integration_glfixed_table_alloc(2*n-1);

  // std::cout << std::chrono::duration_cast<std::chrono::duration<double>>
  //   (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;

  // pre = std::chrono::high_resolution_clock::now();
  // for (int i : boost::irange(0, n)) {
  //   gsl_integration_glfixed(&F, 0, 1, w_table);
  // }
  // std::cout << std::chrono::duration_cast<std::chrono::duration<double>>
  //   (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;

  // pre = std::chrono::high_resolution_clock::now();
  // gsl_integration_workspace * w
  //   = gsl_integration_workspace_alloc (1000);
  // std::cout << std::chrono::duration_cast<std::chrono::duration<double>>
  //   (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;

  // double result, error;
  // double expected = -4.0;
  // double alpha = 1.0;

  // pre = std::chrono::high_resolution_clock::now();

  // for (int i : boost::irange(0, n))
  //   gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
  //                         w, &result, &error);

  // std::cout << std::chrono::duration_cast<std::chrono::duration<double>>
  //   (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;

  // gsl_integration_workspace_free (w);
  // gsl_integration_glfixed_table_free (w_table);
