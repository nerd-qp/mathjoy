/*
  Numerical quadrature

  Relying on gsl

  Two dimension numerical quadrature now shipping with square and trigular
  shape only

  Some keypoint NOTE here is the ability to carry out operation on matrix
 */

#ifndef _QUADRATURE_
#define _QUADRATURE_

/*
  include std libraries
 */
#include <vector>
#include <functional>

// gsl
#include <gsl/gsl_integration.h>

// boost ublas
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

// my own utils lib
#include "utils/utils.hpp"

#define DEFAULT_OROER                     1
#define X_offset                          0
#define Y_offset                          1
#define WEIGHT                            2
#define OFFSET                            3
#define CANONICAL_TRIANGULAR_SPACE        0.5

/*
  For canonical square region
 */
#define CANONICAL_SQUARE_LEFT             -1
#define CANONICAL_SQUARE_RIGHT            1

/*
  For 1 dimension
 */
#define ONE_DIM_OFFSET                    2
#define ONE_DIM_X                         0
#define ONE_DIM_WEIGHT                    1

class NotImplemented : std::logic_error {
public:
  NotImplemented() : std::logic_error("Function not yet implemented!\n") {};
};

namespace one_dimension {

};

namespace two_dimension {

  namespace canonical {
    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::vector;
    class base {
    protected:
      int _degree;

    public:
      // Already defined
      // base() : _degree(DEFAULT_ORDER) {};
      base(int degree = DEFAULT_OROER) : _degree(degree) {};

      // = 0 -> with definition compilation failed
      //virtual const std::vector<double>& get_weight_poly() const = 0;
      virtual const std::vector<double>&
      get_weight_poly() const = 0;

      virtual double sum_up
      (const std::function<double (double, double)> &f) const;

      virtual double
      quadrature_compute
      (const std::function<double (double, double)> &f) const = 0;

      // Operation on matrix
      virtual matrix<double>
      quadrature_compute
      (const matrix<utils::function> &f_matrix) const = 0;

      virtual vector<double>
      quadrature_compute
      (const vector<utils::function> &f_matrix) const = 0;
    };

    class square : public base {
    private:
      // Get the weight and location on [-1, 1]
      // of Gussian quadrature
      static std::vector<double>
      __one_dimension_wx(int degree);

    public:
      const std::vector<double> poly_weight;

      square(int degree = DEFAULT_OROER);

      const std::vector<double>&
      get_weight_poly() const override;

      double quadrature_compute
      (const std::function<double (double, double)> &f) const override;

      matrix<double>
      quadrature_compute
      (const matrix<utils::function> &f_matrix) const override;

      vector<double>
      quadrature_compute
      (const vector<utils::function> &f_matrix) const override;
    };

    class triangular : public base {
    public:

      //triangular() : _degree(DEFAULT_OROER) {};
      triangular(int degree = DEFAULT_OROER) : base(degree) {};

      const std::vector<double>&
      get_weight_poly() const override;

      double quadrature_compute
      (const std::function<double (double, double)> &f) const override;

      matrix<double>
      quadrature_compute
      (const matrix<utils::function> &f_matrix) const override;

      vector<double>
      quadrature_compute
      (const vector<utils::function> &f_matrix) const override;

      const std::vector<std::vector<double>> poly_weight =
        {
         {
          0.33333333333333,    0.33333333333333,    1.00000000000000
         },

         {
          0.16666666666667,    0.16666666666667,    0.33333333333333,
          0.16666666666667,    0.66666666666667,    0.33333333333333,
          0.66666666666667,    0.16666666666667,    0.33333333333333
         },

         {
          0.33333333333333,    0.33333333333333,   -0.56250000000000,
          0.20000000000000,    0.20000000000000,    0.52083333333333,
          0.20000000000000,    0.60000000000000,    0.52083333333333,
          0.60000000000000,    0.20000000000000,    0.52083333333333
         },

         {
          0.44594849091597,    0.44594849091597,    0.22338158967801,
          0.44594849091597,    0.10810301816807,    0.22338158967801,
          0.10810301816807,    0.44594849091597,    0.22338158967801,
          0.09157621350977,    0.09157621350977,    0.10995174365532,
          0.09157621350977,    0.81684757298046,    0.10995174365532,
          0.81684757298046,    0.09157621350977,    0.10995174365532,
         },

         {
          0.33333333333333,    0.33333333333333,    0.22500000000000,
          0.47014206410511,    0.47014206410511,    0.13239415278851,
          0.47014206410511,    0.05971587178977,    0.13239415278851,
          0.05971587178977,    0.47014206410511,    0.13239415278851,
          0.10128650732346,    0.10128650732346,    0.12593918054483,
          0.10128650732346,    0.79742698535309,    0.12593918054483,
          0.79742698535309,    0.10128650732346,    0.12593918054483
         },

         {
          0.24928674517091,    0.24928674517091,    0.11678627572638,
          0.24928674517091,    0.50142650965818,    0.11678627572638,
          0.50142650965818,    0.24928674517091,    0.11678627572638,
          0.06308901449150,    0.06308901449150,    0.05084490637021,
          0.06308901449150,    0.87382197101700,    0.05084490637021,
          0.87382197101700,    0.06308901449150,    0.05084490637021,
          0.31035245103378,    0.63650249912140,    0.08285107561837,
          0.63650249912140,    0.05314504984482,    0.08285107561837,
          0.05314504984482,    0.31035245103378,    0.08285107561837,
          0.63650249912140,    0.31035245103378,    0.08285107561837,
          0.31035245103378,    0.05314504984482,    0.08285107561837,
          0.05314504984482,    0.63650249912140,    0.08285107561837
         }

        };
    };

    // takes a f(x, y) and three coordinates, and transform it into a function
    // that can later be used to compute its quadrature.
    // jacobian and some linear operation
    template <typename FloatT = double> // floating point type
    utils::function function2Canonical(const utils::function& f,
                                       const std::pair<FloatT, FloatT>& p1,
                                       const std::pair<FloatT, FloatT>& p2,
                                       const std::pair<FloatT, FloatT>& p3) {
      FloatT x1, x2, x3, y1, y2, y3;
      std::tie(x1, y1) = p1;
      std::tie(x2, y2) = p2;
      std::tie(x3, y3) = p3;

      // the jacobian function is actually constant
      auto Jacobian = fabs((x2 - x1 ) * (y3 - y1) - (x3 - x1) * (y2 - y1));

      utils::function transformedFun
        {
         //return
         [Jacobian, x1, x2, x3, y1, y2, y3, f]
         (FloatT kexi, FloatT yita) {
           auto x = x1 + (x2 - x1) * kexi + (x3 - x1) * yita;
           auto y = y1 + (y2 - y1) * kexi + (y3 - y1) * yita;

           // std::cout << x << ',' << y << "-> ";
           // std::cout << f(x, y) << '\n';
           return f(x, y) * Jacobian;
         }};
      return transformedFun;
    }
  };

};

#endif
