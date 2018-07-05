#include "quadrature.hpp"
#include <iostream>

namespace two_dimension {
  namespace canonical {
    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::vector;

    // Methods of base
    double base::sum_up
    (const std::function<double (double, double)> &f) const {
      std::vector<double> weight_poly = get_weight_poly();
      double res = 0;

      for ( int i = 0; i < weight_poly.size(); i += OFFSET ) {
        double weight = weight_poly[WEIGHT + i];
        double x = weight_poly[X_offset + i];
        double y = weight_poly[Y_offset + i];

        res += weight * f(x, y);
      }

      return res;
    }

    // Methods of square
    std::vector<double>
    square::__one_dimension_wx(int degree) {
      // TODO
      // This part could be optimized by switching to
      // compute through compilation(pre-compiled compute
      // weight + location)

      // TODO
      // These loops should be destroyed(replaced)
      // And no check bounds here! which is dangerous!!!!

      gsl_integration_glfixed_table * fixed_table =
        gsl_integration_glfixed_table_alloc(degree);

      // result vector and transfer vector
      // temporarily containing weights and nodes
      std::vector<double> res(degree * degree * OFFSET);
      std::vector<double> transfer(degree * ONE_DIM_OFFSET);

      for ( std::size_t i = 0; i < degree; ++ i ) {
        double * x = &transfer[ONE_DIM_OFFSET * i + ONE_DIM_X];
        double * weight = &transfer[ONE_DIM_OFFSET * i + ONE_DIM_WEIGHT];
        gsl_integration_glfixed_point(CANONICAL_SQUARE_LEFT,
                                      CANONICAL_SQUARE_RIGHT,
                                      i, x, weight, fixed_table);
      }

      // destroy table here
      gsl_integration_glfixed_table_free(fixed_table);

      // Generate with tensor product
      for ( int i = 0; i < degree; ++ i ) {
        for ( int j = 0; j < degree; ++ j ) {
          double *x = &res[(i + j * degree) * OFFSET + X_offset];
          double *y = &res[(i + j * degree) * OFFSET + Y_offset];
          double * weight = &res[(i + j * degree) * OFFSET + WEIGHT];

          *x = transfer[ONE_DIM_OFFSET * i + ONE_DIM_X];
          *y = transfer[ONE_DIM_OFFSET * j + ONE_DIM_X];
          *weight =
            transfer[ONE_DIM_OFFSET * j + ONE_DIM_WEIGHT] *
            transfer[ONE_DIM_OFFSET * i + ONE_DIM_WEIGHT];
        }
      }

      return res;
    }

    square::square(int degree) :
      base(degree), poly_weight(__one_dimension_wx(degree)) {}

    const std::vector<double>&
    square::get_weight_poly() const {
      // if ( _degree >= (int)poly_weight.size() )
      //   throw "Degree too high, not supported!\n";
      // else if ( _degree <= 0 )
      //   throw "Degree negative, can't proceed!\n";
      return poly_weight;
    };

    double square::quadrature_compute
    (const std::function<double (double, double)> &f) const {
      return sum_up(f);
    }

    matrix<double>
    square::quadrature_compute
    (const matrix<utils::function> &f_matrix) const {
      int size1 = f_matrix.size1();
      int size2 = f_matrix.size2();

      matrix<double> res(size1, size2);

      for ( auto i = 0; i < size1; ++ i )
        for ( auto j = 0; j < size2; ++ j )
          res(i, j) = quadrature_compute( f_matrix (i, j) );

      return res;
    }

    vector<double>
    square::quadrature_compute
    (const vector<utils::function> &f_vector) const {
      size_t size1 = f_vector.size();

      vector<double> res(size1);

      for ( size_t i = 0; i < size1; ++ i )
        res(i) = quadrature_compute( f_vector(i) );

      return res;
    }

    // Methods of triangle
    const std::vector<double>&
    triangular::get_weight_poly() const {
      if ( _degree >= (int)poly_weight.size() )
        throw "Degree too high, not supported!\n";
      else if ( _degree <= 0 )
        throw "Degree negative, can't proceed!\n";
      return poly_weight[_degree-1];
    };

    double
    triangular::quadrature_compute
    (const std::function<double (double, double)> &f) const {
      return CANONICAL_TRIANGULAR_SPACE * sum_up(f);
    }

    matrix<double>
    triangular::quadrature_compute
    (const matrix<utils::function> &f_matrix) const {
      int size1 = f_matrix.size1();
      int size2 = f_matrix.size2();

      matrix<double> res(size1, size2);

      for ( auto i = 0; i < size1; ++ i )
        for ( auto j = 0; j < size2; ++ j )
          res(i, j) = quadrature_compute( f_matrix (i, j) );

      return res;
    }

    // TODO, these codes can be inherited instead
    vector<double>
    triangular::quadrature_compute
    (const vector<utils::function> &f_vector) const {
      size_t size1 = f_vector.size();

      vector<double> res(size1);

      for ( size_t i = 0; i < size1; ++ i )
        res(i) = quadrature_compute( f_vector(i) );

      return res;
    }

  }
}
