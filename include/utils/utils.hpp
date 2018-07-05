
/*
 * Library for some useful utilities:
 * make_unique, zip
 * std::function wrapped with operator*, +, / overloaded
 */

#ifndef _UTILS_
#define _UTILS_
#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <boost/numeric/ublas/vector.hpp>

namespace utils {
  using boost::numeric::ublas::vector;

  /*
    Implement the almost_equal mostly for floating point comparason
  */
  template<class T>
  typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
  almost_equal(T x, T y, int ulp = 1)
  // as in how many unit of precision, let's default to 1
  {
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
      // unless the result is subnormal
      || std::abs(x-y) < std::numeric_limits<T>::min();
  }


  // TODO could try to extend this
  template<class InputIt, class ForwardIt,
           class op>
  void zip(InputIt first, InputIt last,
           ForwardIt s_first, ForwardIt s_last,
           op f) {
    for (; first != last && s_first != s_last; ++ first, ++ s_first )
      f(*first, *s_first);
  }

  template<typename T, typename... Args>
  std::unique_ptr<T> make_unique(Args&&... args)
  {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  template<typename Type1, typename Type2>
  Type1 inner_prod(const boost::numeric::ublas::vector<Type1> &v1,
                   const boost::numeric::ublas::vector<Type2> &v2) {
    if ( v1.size() == v2.size() && v1.size() > 0) {
      auto base = v1[0] * v2[0];

      for ( size_t i = 1; i < v1.size(); ++ i )
        base = base + v1[i] * v2[i];

      return base;
    }
  }

  /*
   * Term stolen from deep learning
   * It's just operate on all possible combination
   */
  template<typename Vector, typename Functor>
  void fully_connect(const Vector &vec1,
                     const Vector &vec2,
                     const Functor &functor) {

    for ( auto i = vec1.begin(); i != vec1.end(); ++ i )
      for ( auto j = vec2.begin(); j != vec2.end(); ++ j )
        functor(*i, *j);

  }


  // Some function to deal with interpolation with Lagrange
  inline double
  compute_denominator_at(const std::vector<double> &_local_node,
                         size_t index);

  inline std::vector<double>
  compute_denominator_all(const std::vector<double> &_local_node);

  inline std::function<double (double)>
  LagrangePolynomial_at(const std::vector<double> &_local_node,
                     const std::vector<double> &_denominator,
                     size_t index);

  inline std::function<double (double)>
  LagrangePolynomial_diff_at(const std::vector<double> &_local_node,
                          const std::vector<double> &_denominator,
                          int index);

  inline std::vector<std::function<double (double)>>
  LagrangePolynomial(const std::vector<double> &_local_node,
                     const std::vector<double> &_denominator,
                     size_t index);

  inline std::vector<std::function<double (double)>>
  LagrangePolynomial_diff(const std::vector<double> &_local_node,
                          const std::vector<double> &_denominator,
                          int index);

  inline std::vector<double>
  interpolation_gen(double a, double b,
                    const std::vector<double> &hi);

  // inline std::tuple<vector<function>,
  //                   vector<function>,
  //                   vector<function>>
  // Lagrange_with_diff_2dim(const std::vector<double> &x_coor,
  //                         const std::vector<double> &y_coor);
}

#include "utils_function.hpp"

#include "utils_impl.hpp"

#endif
