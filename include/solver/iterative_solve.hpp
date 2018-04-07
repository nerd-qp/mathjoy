#include <boost/numeric/ublas/matrix.hpp>

/*
  A Gauss Seidel Method to solve linear equations
 */

template<typename Matrix, typename Vector>
bool Gauss_Seidel_Method(const Matrix & A,
                         // An initial x should be provided,
                         // no default value
                         // please don't pass in the b==x, it is not allowed!!!
                         const Vector & b, Vector & x,
                         // with low precision like 1e-6 is no good
                         double precision = 1e-20, size_t MAX_ITER = 1e5) {
  // This is an iterative process

  // This variable is used to determine if the process has converged
  double norm = precision + 1;
  size_t count = 0;
  while ( norm > precision && count++ < MAX_ITER) {
    norm = 0;
    for ( size_t i = 0; i < A.size1(); ++ i ) {

      double next_x_i = b(i);

      for ( size_t j = 0; j < i; ++ j )
        next_x_i -= A(i, j) * x(j);

      for ( size_t j = i+1; j < A.size2(); ++ j )
        next_x_i -= A(i, j) * x(j);

      next_x_i /= A(i, i);

      // 1-Norm for absolute error
      // norm += fabs(next_x_i - x(i));

      // infinite-Norm for absolute error
      double error = fabs(next_x_i - x(i));
      norm = norm > error ? norm : error;

      x(i) = next_x_i;
    }
    // if ( count % static_cast<int>(1e4) == 0 )
    //   std::cout << norm << '\n';
  }

  return true;
}
