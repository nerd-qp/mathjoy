#include "utils/matrix_io.hpp"
namespace utils {

  using namespace boost::spirit;
  using qi::double_;
  using qi::phrase_parse;
  using ascii::space;

  // have to say that this is very powerful indeed! In such short lines of code!
  template <typename Iterator>
  bool parse_matrix(Iterator &first, Iterator &last, parseFormat format) {

    bool r = false;
    switch (format)
      {
        // Parse matrix like [1,2;3,4]
      case MATLAB : r = phrase_parse(first,
                                     last,
                                     ('[' >>
                                      ((double_ % ',') % ';')
                                      >> ']'),
                                     space); break;
        // Parse matrix like [1,2;3,4]
      case NUMPY : r = phrase_parse(first,
                                    last,
                                    ('[' >>
                                     ('[' >> double_ % ',' >> ']') %','
                                     >> ']'),
                                    space); break;
      default : {
        return false;
      }
      }

    // doesn't have operator!=, the workaround, using non const iterator
    if (first != last)
      return false;

    return r;
  }

  // for matrix have row, column >= 0, and same columns
  template <typename Iterator> // rightly assumes to be vector<vector>
  inline bool _error_checking_(Iterator &first, Iterator &last, std::vector<std::vector<double>> &matrix) {
    // doesn't have operator!=, the workaround
    if (first != last)
      return false;

    // to check if matrix have same columns
    int row_size = -1;
    if (matrix.size() <= 0 || (row_size = matrix[0].size()) <= 0)
      return false;

    for (const auto &row : matrix)
      if (row_size != row.size())
        return false;

    return true;
  }

  template <typename Iterator>
  bool input_matrix(Iterator &first, Iterator &last, std::vector<std::vector<double>> &matrix, parseFormat format) {
    bool r = false;

    switch (format)
      {
        // Parse matrix like [1,2;3,4]
      case MATLAB : r = phrase_parse(first,
                                     last,
                                     ('[' >>
                                      ((double_ % ',') % ';')
                                      >> ']'),
                                     space,
                                     matrix); break;
        // Parse matrix like [1,2;3,4]
      case NUMPY : r = phrase_parse(first,
                                    last,
                                    ('[' >>
                                     ('[' >> double_ % ',' >> ']') %','
                                     >> ']'),
                                    space,
                                    matrix); break;
      default : {
        return false;
      }
      }

    if (!_error_checking_(first, last, matrix)) {
      return false;
    }

    return r;
  }

  // using matrix = boost::numeric::ublas::matrix<double>;
  // template <typename Iterator>
  // bool input_matrix_numpy(Iterator &first, Iterator &last, matrix& inputMat) {
  //   std::vector<std::vector<double>> vecMat;
  //   if ( input_matrix_numpy(first, last, vecMat) ) {
  //     // because of error check
  //     int row = vecMat.size(), column = vecMat[0].size();
  //     inputMat = matrix(row, column);

  //     for (int i = 0; i < row; ++ i)
  //       for (int j = 0; j < column; ++ j)
  //         inputMat (i, j) = vecMat[i][j];

  //     return true;
  //   }
  //   else return false;
  // }

  using matrix = boost::numeric::ublas::matrix<double>;
  template <typename Iterator>
  bool input_matrix(Iterator &first, Iterator &last, matrix& inputMat, parseFormat format) {
    std::vector<std::vector<double>> vecMat;
    if ( input_matrix(first, last, vecMat, format) ) {
      // because of error check
      int row = vecMat.size(), column = vecMat[0].size();
      inputMat = matrix(row, column);

      for (int i = 0; i < row; ++ i)
        for (int j = 0; j < column; ++ j)
          inputMat (i, j) = vecMat[i][j];

      return true;
    }
    else return false;
  }

  // Initialize templates
  using Iterator = boost::spirit::istream_iterator;
  template bool parse_matrix<Iterator>(Iterator&, Iterator&, parseFormat);
  template bool input_matrix<Iterator>(Iterator &, Iterator &, std::vector<std::vector<double>> &, parseFormat);
  template bool input_matrix<Iterator>(Iterator &, Iterator &, matrix&, parseFormat);

}
