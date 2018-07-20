
#ifndef _MATHJOY_GAUSSIAN_ELIMINATION_HPP_
#define _MATHJOY_GAUSSIAN_ELIMINATION_HPP_

namespace mathjoy {

  /*
    This function takes a matrix and a vector, and use the guassian elimination
    method to solve the linear equation
   */
  template <typename Matrix, typename Vector>
  void Gaussian_elimination(Matrix &A, Vector &b) {
    auto size1 = A.size1();
    auto size2 = A.size2();

    for ( auto pivot_pos = 0ul; pivot_pos < size1; ++ pivot_pos ) {

      // NOTE float? anyway to make it more flexible?
      auto pivot = static_cast<double>(A(pivot_pos, pivot_pos));
      // Stop when pivot is zero
      if ( pivot == 0 ) return;

      for ( auto elim_row = pivot_pos + 1; elim_row < size1; ++ elim_row ) {

        auto first_elem = A(elim_row, pivot_pos);
        // no operation needed if already 0
        if ( first_elem == 0 ) break;

        auto m = first_elem / pivot;
        // eliminate the first column of the row
        A(elim_row, pivot_pos) = 0;
        b(elim_row) -= b(pivot_pos) * m;

        for ( auto col = pivot_pos + 1; col < size2; ++ col )
          A(elim_row, col) -= A(pivot_pos, col) * m;

      }

    }
  }

}

#endif
