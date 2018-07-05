
/*
  Provides the basic input output methods for custom input file
 */

#ifndef _MATRIX_IO_
#define _MATRIX_IO_
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <vector>

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

#endif
