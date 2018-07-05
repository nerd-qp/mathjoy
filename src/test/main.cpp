

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
#define MAXIMUM_NODES 150

using boost::numeric::ublas::matrix;
using boost::numeric::ublas::vector;

void for_testing(std::string message,
                 const std::vector<double> &xi,
                 const std::vector<double> &hi,
                 const matrix<double> &stiffness_matrix,
                 const vector<double> &right_f_phi,
                 const vector<double> &solution,
                 int degree = 1)
{

  using std::endl;
  using std::cout;

  cout << message << endl;

  // calculate error here
  // cout << stiffness_matrix << endl

  //      << endl

  //      << right_f_phi << endl

  //      << endl

  //      << solution << endl;

  double total_error = 0.;
  double _max = -1;

  // dealing with index 0, treating as a special case
  std::vector<double> local_x_i(degree);
  for ( const auto &i : boost::irange(1, degree) ) {
    local_x_i[i] = i * hi[0] / (double)degree + xi[0];
    total_error += fabs(result_f(local_x_i[i]) - solution(i-1));
    _max = std::max(fabs(result_f(local_x_i[i]) - solution(i-1)), _max);
  }

  for ( const auto &i : boost::irange(1, (int)xi.size()-1) ) {

    for ( const auto &j : boost::irange(0, degree) )
      local_x_i[j] =
        xi[i] + j * hi[i] / (double)degree;

    for ( const auto &j : boost::irange(0, degree) ) {

      _max = std::max(fabs(result_f(local_x_i[j]) -
                           solution(i * degree + j - 1)), _max);
      total_error +=
        fabs(result_f(local_x_i[j]) - solution(i * degree + j - 1));

      //std::cout << result_f(local_x_i[j]) << "\n";
      //std::cout << i * degree + j - 1 << "\n";

    }
    //cout << total_error << endl;
  }

  // dealing with index n-1, treating as a special case
  _max =
    std::max(fabs(result_f(xi[xi.size()-1]) -
             solution(solution.size()-1)), _max);

  total_error += fabs(result_f(xi[xi.size()-1]) -
                      solution(solution.size() - 1));

  std::cout << message << endl
            << "Interpolation nodes: \n"
            << solution.size()
            << "\n"

            << "Total error is : \n"
            << total_error
            << "\n"

            << "Mean error is : \n"
            << total_error / (double)solution.size()
            << "\n"

            << "Max error is : \n"
            << _max
            << "\n"

            // << "Total error divided by hi square is : \n"
            // << total_error / (double)(solution.size() * solution.size())
            // << "\n"
            << endl;
}

int main() {

  for ( int i = 4; i < MAXIMUM_NODES; i *= 2 ) {
    std::cout << "starting general\n\n";
    // General solver
    using namespace general;
    //int i = 1;
    int b = 1; int a = 0;

    double degree = 2;

    {

      std::vector<double> hi(i/degree, degree*(b-a)/(double)i), xi;
      FEM_solver::x_i_array(a, b, hi, xi);

      auto pre = std::chrono::high_resolution_clock::now();
      FEM_solver solver_general(a, b, hi, &test_p, &test_q, &test_f, degree);
      std::cout << "Time used to initialize: " <<
        std::chrono::duration_cast<std::chrono::duration<double>>
        (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;

      auto stiffness_matrix = solver_general.get_stiffness_matrix();
      auto right_f_phi = solver_general.get_right_f_phi();
      auto solution = solver_general.get_solution();

      std::cout << "Total time used: " <<
        std::chrono::duration_cast<std::chrono::duration<double>>
        (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;

      for_testing("General",
                  xi, hi,
                  stiffness_matrix,
                  right_f_phi,
                  solution,
                  degree);

    }

    {
      std::cout << "starting linear solver\n\n";

      // Linear solver
      using namespace linear;
      std::vector<double> hi(i, (b-a)/(double)i), xi;
      finite_element_solver::x_i_array(a, b, hi, xi);
      auto pre = std::chrono::high_resolution_clock::now();
      finite_element_solver solver_linear(a, b, hi, &test_p, &test_q, &test_f);
      std::cout << "Time used to initialize: " <<
        std::chrono::duration_cast<std::chrono::duration<double>>
        (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;


      auto stiffness_matrix = solver_linear.get_stiffness_matrix();
      auto right_f_phi = solver_linear.get_right_f_phi();
      auto solution = solver_linear.get_solution();

      std::cout << "Total time used: " <<
        std::chrono::duration_cast<std::chrono::duration<double>>
        (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;

      for_testing("Linear",
                  xi, hi,
                  stiffness_matrix,
                  right_f_phi,
                  solution,
                  1);

    }

  }

}
