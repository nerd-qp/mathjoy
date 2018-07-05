#include "quadrature.hpp"
#include "mesh_generator/mesh_generator.hpp"
#include "FEM_2dim.hpp"
#include "utils/utils.hpp"

#include <cmath>
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


using boost::numeric::ublas::matrix;
using boost::numeric::ublas::vector;

void test_square(const utils::function &p,
                 const utils::function &q,
                 const utils::function &f,
                 const utils::function &target_f,
                 size_t maximum) {

  for ( auto i = 1; i <= maximum; i = i << 1 ) {
    linear::FEM_square test_again(std::make_pair<double, double>(-1, 1),
                                  std::make_pair<double, double>(-1, 1),
                                  i, i, p, q, f, 1);

    // std::cout << test_again.get_stiffness_matrix() << "\n\n";
    // std::cout << test_again.get_right_f_phi() << std::endl;

    auto tmp = test_again.get_solution();
    std::cout << tmp(tmp.size() / 2) << "\n";

    double max_diff = 0;
    size_t max_index = -1;
    double sum = 0;

    auto vertices = test_again.get_interpolation_vertices();

    for ( size_t i = 0; i < vertices.size(); ++ i ) {

      auto results = target_f(vertices[i].getX(),
                              vertices[i].getY());

      // std::cout << vertices[i].getX() << ", "
      //           << vertices[i].getY() << " -> "
      //           << tmp(i) << " vs "
      //           << results << '\n';

      auto diff = fabs(results - tmp(i));

      sum += max_diff;

      if ( diff > max_diff ) {
        max_index = i;
        max_diff = diff;
      }

      // if ( min_diff > diff ) {
      //   min_diff = diff;
      // }

    }

    std::cout << "Maximum difference is: "
              << max_diff// << ", "
      //<< vertices[max_index].getX() << ", "
      //<< vertices[max_index].getY() << ", "
      //      << max_index
              << '\n';
    std::cout << "Average difference is: "
              << sum / (double)vertices.size() << '\n';
    // std::cout << "Minimum difference is: "
    //           << min_diff << '\n';
  }

}

void test_triangular(const utils::function &p,
                     const utils::function &q,
                     const utils::function &f,
                     const utils::function &target_f,
                     size_t maximum) {
  for ( auto i = 1; i <= maximum; i = i << 1 ) {
    linear::FEM_triangular test_again(std::make_pair<double, double>(-1, 1),
                                      std::make_pair<double, double>(-1, 1),
                                      i, i, p, q, f, 1);

    // std::cout << test_again.get_stiffness_matrix() << "\n\n";
    // std::cout << test_again.get_right_f_phi() << std::endl;

    auto tmp = test_again.get_solution();
    std::cout << tmp(tmp.size() / 2) << "\n";

    double max_diff = 0;
    size_t max_index = -1;
    double sum = 0;

    auto vertices = test_again.get_interpolation_vertices();

    for ( size_t i = 0; i < vertices.size(); ++ i ) {

      auto results = target_f(vertices[i].getX(),
                              vertices[i].getY());

      // std::cout << vertices[i].getX() << ", "
      //           << vertices[i].getY() << " -> "
      //           << tmp(i) << " vs "
      //           << results << '\n';

      auto diff = fabs(results - tmp(i));

      sum += max_diff;

      if ( diff > max_diff ) {
        max_index = i;
        max_diff = diff;
      }

      // if ( min_diff > diff ) {
      //   min_diff = diff;
      // }

    }

    std::cout << "Maximum difference is: "
              << max_diff// << ", "
      //<< vertices[max_index].getX() << ", "
      //<< vertices[max_index].getY() << ", "
      //      << max_index
              << '\n';
    std::cout << "Average difference is: "
              << sum / (double)vertices.size() << '\n';
    // std::cout << "Minimum difference is: "
    //           << min_diff << '\n';
  }

}

int main() {
  using namespace two_dimension::canonical;

  utils::function p([=] (double x, double y) -> double {
                      return 1;
                      //return 0;
                    });

  utils::function q([=] (double x, double y) -> double {
                      //return 1;
                      return 0;
                     });

  utils::function f([=] (double x, double y) -> double {
                      //return (x*x - 1) * (y*y-1);
                      //return x*x - 1;
                      //return 1;
                      //return M_PI * M_PI * sin(M_PI * x);
                      //return -2;
                      return - (2 * y * y + 2 * x * x - 4);
                      // return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
                      // return sin(M_PI * x) * sin(M_PI * y);
                    });

  utils::function target_f([=] (double x, double y) -> double {
                             //return 1;
                             //return sin(M_PI * x);
                             //return x * x - 1;
                             //return 0;
                             return (x*x - 1) * (y*y - 1);
                             //return sin(M_PI * x) * sin(M_PI * y);
                           });

  test_square(p, q, f, target_f, 32);

  std::cout << "Next up is triangular shape FEM" << std::endl;

  test_triangular(p, q, f, target_f, 32);

}
