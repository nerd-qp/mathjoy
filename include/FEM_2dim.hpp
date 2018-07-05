/*
  The core of 2 dimensional FEM solver

  stiffness matrix assembly
  stiffness matrix solving

*/

/*
  include self-made libraries
*/
#include "mesh_generator/mesh_generator.hpp"
#include "quadrature.hpp"
#include "utils/utils.hpp"

/*
  Boost matrix libraries
*/
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

/*
  Include std libraries
*/
#include <memory>

#ifndef _FINITE_ELEMENT_METHOD_TWO_DIMENSION_
#define _FINITE_ELEMENT_METHOD_TWO_DIMENSION_

/*
 * Some configuration to enable move semantics
 */
#define BOOST_UBLAS_MOVE_SEMANTICS

namespace linear {

  using boost::numeric::ublas::matrix;
  using boost::numeric::ublas::vector;
  using boost::numeric::ublas::mapped_matrix;

  // Basics
  class FEM_solver_2dim_base {
  protected:
    mapped_matrix<double>  _stiffness_matrix;
    vector<double>         _right_f_phi;
    std::unique_ptr<simple_mesh> _mesh;

    utils::function _p;
    utils::function _q;
    utils::function _f;

  public:
    // Give 3 or 4 vertex(a face) returns
    // a matrix of (function, index) pair to compute
    // index for assembling
    virtual std::pair<matrix<utils::function>,
                      matrix<std::pair<size_t, size_t>>>
    K_element_matrix(const face &element_face) const = 0;

    virtual std::pair<matrix<utils::function>,
                      matrix<std::pair<size_t, size_t>>>
    M_element_matrix(const face &element_face) const = 0;

    virtual std::pair<vector<utils::function>, vector<size_t>>
    f_element_vector(const face &element_face) const = 0;

    /*
     * Interpolation funciton generation
     */
    // inline vector<utils::function>
    // Lagrange_gen(size_t degree) const;

    // // this thing is buggy... return value is wrong
    // // some sort of copy error
    // inline std::pair<vector<utils::function>,
    //                  vector<utils::function>>
    // Lagrange_diff_gen(size_t degree) const;

    /*
     * assemble stiffness matrix and right function vector
     * Since the functions generated from the above functions
     * are just temporary values, therefore no need to store
     * them
     */
    virtual void
    assemble_stiffness_matrix() = 0;

    virtual void
    assemble_right_f_phi() = 0;

    virtual void
    handle_boundary_point(size_t index) = 0;

  public:

    // don't know how to express region
    // just going to provide p, q, f initialization
    // FEM equation N is determined as well
    FEM_solver_2dim_base(utils::function p,
                         utils::function q,
                         utils::function f,
                         size_t size_N) :
      _stiffness_matrix(size_N, size_N),
      _right_f_phi(size_N, 0), _p(p), _q(q), _f(f)
    {  };

    inline mapped_matrix<double>
    get_stiffness_matrix() const {
      return _stiffness_matrix;
    }

    inline vector<double>
    get_right_f_phi() const {
      return _right_f_phi;
    };

    inline const std::vector<vertex> &
    get_interpolation_vertices() const {
      return _mesh->get_vertices();
    };

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

  };

  /*
   * This class only deals with the simple square region
   * Require manual region specification
   */
  class FEM_square : public FEM_solver_2dim_base {
  private:
    // The base function vector
    vector<utils::function> _base_fun_vec;
    // The base function differentiation vector
    vector<utils::function> _base_diff_vec_kexi;
    vector<utils::function> _base_diff_vec_yita;
    std::vector<double>     _denominator;
    std::unique_ptr<two_dimension::canonical::square> _quadrature_table;

    // void
    // Lagrange_diff_gen(size_t degree,
    //                   vector<utils::function> *diff_kexi,
    //                   vector<utils::function> *diff_yita);

    /*
     * A linear coordinate transformation
     * Take two vector of x, y coordinates
     * return a pair of function that maps kexi, yita
     * from [-1, 1] x [-1, 1] to the corresponding region
     * represented by x, y coordinates
     */
    inline std::pair<utils::function, utils::function>
    __linear_co_transform(const vector<double> &x_coor,
                          const vector<double> &y_coor) const;

    /*
     * Compute the derivative function for x, y coordinates
     */
    inline std::tuple<utils::function,
                      utils::function,
                      utils::function,
                      utils::function>
    __local_derivative_function_get(const vector<double> &x_coor,
                                    const vector<double> &y_coor) const;

    /*
     * Some helper funciton to get the corresponding
     * x_coor, y_coor and vertices_index
     * since multiple places involved
     * Please forgive me for the name...
     * personal use only...
     */

    inline std::tuple<std::vector<size_t>, vector<double>, vector<double>>
    __index_x_coor_y_coor_get(const face &element_face) const;

    /*
     * As suggested by the name
     * Pre: Have to prepare the mesh and have generated all
     * required base functions
     */
    inline void
    __assemble_local_stiffness_matrix(const face &element_face);

    inline void
    __assemble_local_right_f_phi(const face &element_face);

    /*
     * delegating constructor
     */

    FEM_square(std::pair<double, double> v1,
               std::pair<double, double> v2,
               size_t slice_x_num,
               size_t slice_y_num,
               utils::function p,
               utils::function q,
               utils::function f,
               size_t degree,
               const size_t &freedom);

    /*
     * Pre: mesh generation complete
     */
    void
    assemble_stiffness_matrix() override;

    void
    assemble_right_f_phi() override;

    /*
     * Should only be used post stiffness matrix assembling
     * and right_f_phi assembing done
     */
    inline void
    handle_boundary_point(size_t index) override;

  public:
    FEM_square(std::pair<double, double> v1,
               std::pair<double, double> v2,
               size_t slice_x_num,
               size_t slice_y_num,
               utils::function p,
               utils::function q,
               utils::function f,
               size_t degree = 1);

    std::pair<matrix<utils::function>, matrix<std::pair<size_t, size_t>>>
    K_element_matrix(const face &element_face) const override;

    std::pair<matrix<utils::function>, matrix<std::pair<size_t, size_t>>>
    M_element_matrix(const face &element_face) const override;

    std::pair<vector<utils::function>, vector<size_t>>
    f_element_vector(const face &element_face) const override;

  };

  /*
   * This class only deals with the simple square region
   * Require manual region specification
   */
  class FEM_triangular : public FEM_solver_2dim_base {
  private:
    // The base function vector
    vector<utils::function> _base_fun_vec;
    // The base function differentiation vector
    vector<utils::function> _base_diff_vec_kexi;
    vector<utils::function> _base_diff_vec_yita;
    std::vector<double>     _denominator;
    std::unique_ptr<two_dimension::canonical::triangular> _quadrature_table;

    /*
     * A linear coordinate transformation
     * Take two vector of x, y coordinates
     * return a pair of function that maps kexi, yita
     * from [-1, 1] x [-1, 1] to the corresponding region
     * represented by x, y coordinates
     */
    inline std::pair<utils::function, utils::function>
    __linear_co_transform(const vector<double> &x_coor,
                          const vector<double> &y_coor) const;

    /*
     * Compute the derivative function for x, y coordinates
     */
    inline std::tuple<utils::function,
                      utils::function,
                      utils::function,
                      utils::function>
    __local_derivative_function_get(const vector<double> &x_coor,
                                    const vector<double> &y_coor) const;

    /*
     * Some helper funciton to get the corresponding
     * x_coor, y_coor and vertices_index
     * since multiple places involved
     * Please forgive me for the name...
     * personal use only...
     */

    inline std::tuple<std::vector<size_t>, vector<double>, vector<double>>
    __index_x_coor_y_coor_get(const face &element_face) const;

    /*
     * As suggested by the name
     * Pre: Have to prepare the mesh and have generated all
     * required base functions
     */
    inline void
    __assemble_local_stiffness_matrix(const face &element_face);

    inline void
    __assemble_local_right_f_phi(const face &element_face);

    /*
     * delegating constructor
     */

    FEM_triangular(std::pair<double, double> v1,
                   std::pair<double, double> v2,
                   size_t slice_x_num,
                   size_t slice_y_num,
                   utils::function p,
                   utils::function q,
                   utils::function f,
                   size_t degree,
                   const size_t &freedom);

    /*
     * Pre: mesh generation complete
     */
    void
    assemble_stiffness_matrix() override;

    void
    assemble_right_f_phi() override;

    /*
     * Should only be used post stiffness matrix assembling
     * and right_f_phi assembing done
     */
    inline void
    handle_boundary_point(size_t index) override;

  public:
    FEM_triangular(std::pair<double, double> v1,
                   std::pair<double, double> v2,
                   size_t slice_x_num,
                   size_t slice_y_num,
                   utils::function p,
                   utils::function q,
                   utils::function f,
                   size_t degree = 1);

    std::pair<matrix<utils::function>, matrix<std::pair<size_t, size_t>>>
    K_element_matrix(const face &element_face) const override;

    std::pair<matrix<utils::function>, matrix<std::pair<size_t, size_t>>>
    M_element_matrix(const face &element_face) const override;

    std::pair<vector<utils::function>, vector<size_t>>
    f_element_vector(const face &element_face) const override;

  };


};

#endif
