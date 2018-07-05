/*
 * Intended to generate a LagrangePolynomial on two dimension
 * TODO Or even higher dimension
 */


#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>

namespace utils {

  inline double
  compute_denominator_at(const std::vector<double> &_local_node,
                         size_t index) {
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

  inline std::vector<double>
  compute_denominator_all(const std::vector<double> &_local_node) {
    std::vector<double> denominators(_local_node.size());

    for ( size_t i = 0; i < denominators.size(); ++ i )
      denominators[i] = compute_denominator_at(_local_node, i);

    return denominators;
  }

  /*
   * Generate a Lagrange Polynomial given a intepolation vector
   * and the denominator vector
   * As a helper function
   */
  inline std::function<double (double)>
  LagrangePolynomial_at(const std::vector<double> &_local_node,
                        const std::vector<double> &_denominator,
                        size_t index) {
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

  /*
   *
   */
  inline std::function<double (double)>
  LagrangePolynomial_diff_at(const std::vector<double> &_local_node,
                             const std::vector<double> &_denominator,
                             int index) {
    auto start = _local_node.cbegin();
    auto mid = _local_node.cbegin() + index;
    auto end = _local_node.cend();

    // First nasty version
    return
      ([=] ( double x ) {
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
  }

  /*
   * Given n nodes on some region
   * Generate LagrangePolynomial function
   * at every position of nodes
   */
  inline std::vector<std::function<double (double)>>
  LagrangePolynomial(const std::vector<double> & _local_node) {
    std::vector<std::function<double (double)>> res(_local_node.size());
    std::vector<double> denominators = compute_denominator_all(_local_node);

    for ( size_t i = 0; i < res.size(); ++ i )
      res[i] = LagrangePolynomial_at(_local_node, denominators, i);

    return res;
  }

  /*
   * Given n nodes on some region
   * Generate LagrangePolynomial differentiation function
   * at every position of nodes
   */
  inline std::vector<std::function<double (double)>>
  LagrangePolynomial_diff(const std::vector<double> &_local_node) {
    std::vector<std::function<double (double)>> res(_local_node.size());
    std::vector<double> denominators = compute_denominator_all(_local_node);

    for ( size_t i = 0; i < res.size(); ++ i )
      res[i] = LagrangePolynomial_diff_at(_local_node, denominators, i);

    return res;
  }

  /*
   * Generate intepolation based on a vector of distance h
   * Would throw exception winhen the distance doesn't fully
   * fill the whole region
   */
  inline std::vector<double>
  interpolation_gen(double a, double b,
                    const std::vector<double> &hi)
  {
    double tmp = a;
    std::vector<double> res;
    res.reserve(hi.size()+1);
    for (const auto &h : hi) {
      if ( tmp <= b )
        res.emplace_back(tmp);
      else
        throw "Interpolation out of bound!\n";
      tmp += h;
    }
    res.emplace_back(tmp);

    if ( !almost_equal(tmp, b) )
      throw "Last node is not equal to b!\n";

    return res;
  }

  /*
   * Takes in a set of x, and y coor
   */
  static inline std::tuple<vector<function>,
                           vector<function>,
                           vector<function>>
  Lagrange_with_diff_2dim(const std::vector<double> &x_coor,
                          const std::vector<double> &y_coor) {

    // generate LagrangePolynomial function and
    // differentiation function respectively
    auto tmp_x = LagrangePolynomial(x_coor);
    auto tmp_y = LagrangePolynomial(y_coor);

    auto tmp_diff_x = LagrangePolynomial_diff(x_coor);
    auto tmp_diff_y = LagrangePolynomial_diff(y_coor);

    // copy std::vector -> ublas::vector
    vector<utils::function_1dim> x_onedim_functions(tmp_x.size());
    vector<utils::function_1dim> y_onedim_functions(tmp_y.size());

    vector<utils::function_1dim> x_onedim_diff(tmp_x.size());
    vector<utils::function_1dim> y_onedim_diff(tmp_y.size());

    for ( size_t i = 0; i < x_onedim_functions.size(); ++ i ) {
      x_onedim_functions(i) = tmp_x[i];
      x_onedim_diff(i) = tmp_diff_x[i];
    }

    for ( size_t i = 0; i < y_onedim_functions.size(); ++ i ) {
      y_onedim_functions(i) = tmp_y[i];
      y_onedim_diff(i) = tmp_diff_y[i];
    }

    // not working code

    // std::copy(tmp_x.begin(), tmp_x.end(), x_onedim_functions.begin());
    // std::copy(tmp_y.begin(), tmp_y.end(), y_onedim_functions.begin());

    // std::copy(tmp_diff_x.begin(), tmp_diff_x.end(),
    //           x_onedim_functions.begin());
    // std::copy(tmp_diff_y.begin(), tmp_diff_y.end(),
    //           y_onedim_functions.begin());

    // make two dimensional function with a tensor product
    vector<function> fun_vec =
      tensor_product_fun(x_onedim_functions, y_onedim_functions);

    vector<function> diff_vec_x =
      tensor_product_fun(x_onedim_diff, y_onedim_functions);

    vector<function> diff_vec_y =
      tensor_product_fun(x_onedim_functions, y_onedim_diff);

    // for ( const auto &e : fun_vec ) {
    //   std::cout << e(1, 1) << std::endl;
    // }
    // std::cout << fun_vec[0](10, 10) << std::endl;

    // return
    return std::make_tuple(std::move(fun_vec),
                           std::move(diff_vec_x),
                           std::move(diff_vec_y));
  }
}


  // static inline std::tuple<vector<function>,
//                            vector<function>,
//                            vector<function>>
//   Lagrange_with_diff_2dim(const std::vector<double> &x_coor,
//                           const std::vector<double> &y_coor) {

//     // generate LagrangePolynomial function and
//     // differentiation function respectively
//     auto tmp_x = LagrangePolynomial(x_coor);
//     auto tmp_y = LagrangePolynomial(y_coor);

//     auto tmp_diff_x = LagrangePolynomial_diff(x_coor);
//     auto tmp_diff_y = LagrangePolynomial_diff(y_coor);

//     // copy std::vector -> ublas::vector
//     vector<utils::function_1dim> x_onedim_functions(tmp_x.size());
//     vector<utils::function_1dim> y_onedim_functions(tmp_y.size());

//     vector<utils::function_1dim> x_onedim_diff(tmp_x.size());
//     vector<utils::function_1dim> y_onedim_diff(tmp_y.size());

//     for ( size_t i = 0; i < x_onedim_functions.size(); ++ i ) {
//       x_onedim_functions(i) = tmp_x[i];
//       x_onedim_diff(i) = tmp_diff_x[i];
//     }

//     for ( size_t i = 0; i < y_onedim_functions.size(); ++ i ) {
//       y_onedim_functions(i) = tmp_y[i];
//       y_onedim_diff(i) = tmp_diff_y[i];
//     }

//     // not working code

//     // std::copy(tmp_x.begin(), tmp_x.end(), x_onedim_functions.begin());
//     // std::copy(tmp_y.begin(), tmp_y.end(), y_onedim_functions.begin());

//     // std::copy(tmp_diff_x.begin(), tmp_diff_x.end(),
//     //           x_onedim_functions.begin());
//     // std::copy(tmp_diff_y.begin(), tmp_diff_y.end(),
//     //           y_onedim_functions.begin());

//     // make two dimensional function with a tensor product
//     vector<function> fun_vec =
//       tensor_product_fun(x_onedim_functions, y_onedim_functions);

//     vector<function> diff_vec_x =
//       tensor_product_fun(x_onedim_diff, y_onedim_functions);

//     vector<function> diff_vec_y =
//       tensor_product_fun(x_onedim_functions, y_onedim_diff);

//     // for ( const auto &e : fun_vec ) {
//     //   std::cout << e(1, 1) << std::endl;
//     // }
//     // std::cout << fun_vec[0](10, 10) << std::endl;

//     // return
//     return std::make_tuple(std::move(fun_vec),
//                            std::move(diff_vec_x),
//                            std::move(diff_vec_y));
//   }
// }
