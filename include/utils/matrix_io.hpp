
/*
  Provides the basic input output methods for custom input file
 */

#ifndef _MATHJOY_MATRIX_IO_HPP_
#define _MATHJOY_MATRIX_IO_HPP_

// boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>

// stl
#include <algorithm>
#include <vector>

namespace utils {
  enum parseFormat{MATLAB, NUMPY};

  template <typename Iterator>
  bool parse_matrix(Iterator &first, Iterator &last, parseFormat format);

  template <typename Iterator>
  bool input_matrix(Iterator &, Iterator &,
                    std::vector<std::vector<double>> &,
                    parseFormat);

  template <typename Iterator>
  bool input_matrix(Iterator &, Iterator &,
                    boost::numeric::ublas::matrix<double> &,
                    parseFormat);


  /*
    Accepts input and put into, providing boost::spirit::istream_iterator
    for now vector<vector<double>>

    The vector<vector> passed in should be empty, otherwise, the item would get appended
  */
  // template <typename Iterator, typename Mat_Type>
  // bool input_matrix_numpy(Iterator &first, Iterator &last, Mat_Type &);

  // template <typename Iterator, typename Mat_Type>
  // bool input_matrix_matlab(Iterator &first, Iterator &last, Mat_Type &);

  // std::istream &operator>>(std::istream &_input_stream, matrix<double> & _recv_matrix) {
  //   // use vector to contain


  //   // fail to read the format
  //   _input_stream.setstate(std::ios_base::failbit);
  //   return _input_stream;
  // }

  /*
    if the matrix size and vector(used to initilize the matrix) size is
    inconsistent. It would return a 1 by 1 matrix to indicate failure
  */
  template <typename T, typename F=boost::numeric::ublas::row_major>
  boost::numeric::ublas::matrix<T> makeMatrix(std::size_t m, std::size_t n,
                                              const std::vector<T> &v) {
    using namespace boost::numeric::ublas;
    if ( m * n != v.size() ) return matrix<T>(1, 1, 0);

    unbounded_array<T> storage(m*n);
    std::copy(v.begin(), v.end(), storage.begin());
    return matrix<T>(m, n, storage);
  }


}
#endif
