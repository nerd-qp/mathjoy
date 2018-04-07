#include <boost/numeric/ublas/matrix.hpp>

/*
  A Gauss Seidel Method to solve linear equations
 */

template<typename Matrix, typename Vector, typename number>
bool Gauss_Seidel_Method(const Matrix & A,
                         // An initial x should be provided,
                         // no default value
                         // please don't pass in the b==x, it is not allowed!!!
                         const Vector & b, Vector & x,
                         number precision, size_t MAX_ITER = 1e2) {
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

      norm += fabs(next_x_i - x(i));
      x(i) = next_x_i;
    }
  }

  return true;
}
