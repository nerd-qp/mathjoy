
#ifndef _FVM_2DIM_HPP_
#define _FVM_2DIM_HPP_
/*
  This implementation relies on the external libraries for generating trianglation
 */

/*
  std libraries type traits, before delaunary, because it relies on it...
  some bug that I some fix TODO
*/
#include <type_traits>

/*
  include delaunary generation libraries
 */
#include "delaunay-triangulation/vector2.h"
#include "delaunay-triangulation/triangle.h"
#include "delaunay-triangulation/delaunay.h"

/*
  Boost matrix libraries
 */
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

/*
  need the support of quadrature and utils
 */
#include "utils/utils.hpp"
#include "quadrature.hpp"

namespace FVM {
  using boost::numeric::ublas::vector;
  using boost::numeric::ublas::mapped_matrix;
  class FVM_solver_2dim {

  private:
    mapped_matrix<double>        _stiffness_matrix;
    vector<double>               _right_f_phi;
    std::vector<Vector2<double>> _points;
    // bool mask to determine which point is on boundary
    std::vector<bool>            _on_boundary;
    Delaunay<double>             _triangulation;
    utils::function              _f;
    std::unique_ptr<two_dimension::canonical::triangular> _quadrature_table;

    /*
      compute connectivities of nodes, given triangulation and points
     */
    // void _connect(std::vector<std::vector<std::pair<std::size_t, std::size_t>>>& connect_mapping);
    void _connect(std::vector<std::vector<std::size_t>>& triangles_indices,
                  std::vector<std::vector<std::size_t>>& self_indices,
                  std::vector<std::vector
                  <std::pair<std::size_t, std::size_t>>>& other_indices);

    /*
      Compute all the circum center
     */
    void _circumCenter
    (std::vector<std::vector<std::size_t>>& triangles_indices,
     std::vector<std::vector<std::size_t>>& self_indices,
     std::vector<std::vector
     <std::pair<std::size_t, std::size_t>>>& other_indices,
     std::vector<std::vector<Vector2<double>>>& circumCenter_points);

    /*
      Mark other indices
     */
    // void _markOtherIndices(const std::vector<std::size_t>& self_index,
    //                        <std::pair<std::size_t, std::size_t>>& other_indices);

    /*
      Rearrange triangles given their global indices
     */
    void _reArrange
    (std::vector<std::size_t>& triangles,
     std::vector<std::size_t>& self_indices,
     std::vector<std::pair<std::size_t, std::size_t>>& other_indices);

  public:
    FVM_solver_2dim(const std::vector<Vector2<double>>& points,
                    const Delaunay<double>& triangulation,
                    const std::vector<bool>& on_boundary,
                    const utils::function& f)
      : _points(points), _triangulation(triangulation), _on_boundary(on_boundary),
        _stiffness_matrix(points.size(), points.size()),
        _right_f_phi(points.size()),
        _quadrature_table(utils::make_unique
                          <two_dimension::canonical::triangular>(5)),
        _f(f)
        // assume total nodes N in triangulation
        // the stiffness matrix should be the size NxN
    { };

    /*
      assemble stiffness matrix and right_f_phi
      The algorithm is straightforward

      For every node in the generated triangulation
      find the triangle t that owns it, find in t the other two nodes, and add
      them to set S.

      Compute circumcenter on each triangle

      Compute the weight based on the formula and add to stiffness matrix

      Compute the right side

      done

      return the edges linked by circumCenters of each neighbour triangle of
      nodes
     */
    std::vector<Edge<double>>
    assemble_all();

    inline mapped_matrix<double>
    get_stiffness_matrix() const {
      return _stiffness_matrix;
    }

    inline vector<double>
    get_right_f_phi() const {
      return _right_f_phi;
    }

    inline vector<double>
    get_solution() const {
      auto my_matrix = this->get_stiffness_matrix();
      auto my_vector = this->get_right_f_phi();

      boost::numeric::ublas::permutation_matrix<size_t> pm(my_matrix.size1());
      lu_factorize(my_matrix, pm);
      lu_substitute(my_matrix, pm, my_vector);

      //std::cout << my_vector << std::endl;
      return my_vector;
    }

    inline
    const std::vector<Vector2<double>>& getPoints() {
      return _points;
    }

  };
}
#endif
