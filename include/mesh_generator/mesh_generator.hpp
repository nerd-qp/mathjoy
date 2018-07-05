/*
  This class generates the a simple rectangle mesh area
  and store in a vector-like data structure
  (face -> edge -> vertex)
*/

#ifndef _MESH_GENERATOR_
#define _MESH_GENERATOR_

/*
  std libraries
*/
#include <vector>
#include <memory>
#include "vertex.hpp"
#include "edge.hpp"
#include "face.hpp"
//#include ""

// class mesh {
// protected:

// public:
//   // std::vector<vertex> _global_vertices;
//   // std::vector<edge> _global_edges;
//   // std::vector<face> _global_faces;

// };

class simple_mesh {
protected:
  typedef std::size_t size_t;

  // Although here is a simple rectangle vertex using a
  // matrix representation is better
  // but considering other shape, using vector here

  // NOTE
  // And most importantly using matrix is hard for construction
  // of global stiffness matrix
  std::vector<vertex> _global_vertices;
  std::vector<edge> _global_edges;
  std::vector<face> _global_faces;

public:
  // does nothing
  simple_mesh() {  };

  simple_mesh(size_t vertex_size,
              size_t edge_size,
              size_t face_size) :
    _global_vertices(vertex_size),
    _global_edges(edge_size),
    _global_faces(face_size)
  { };

  const std::vector<face> &
  get_faces() { return _global_faces; };

  const std::vector<vertex> &
  get_vertices() { return _global_vertices; };


};

/*
  Simple mesh only takes care of a simple rectangle region
  from (x1, y1) -> (x2, y2)
*/
class simple_square_mesh : public simple_mesh {
private:
  typedef std::size_t size_t;

  std::size_t _x_num;
  std::size_t _y_num;

  /*
   * locate edge by vertex
   */

  // Link single edges
  inline void
  __link_edges(size_t v1_index,
               size_t v2_index,
               size_t edge_index);

  // link edges recursively

  // NOTE virtual function should not be called in constructors
  // Very dangerous move, see reference below
  /*
   * https://wiki.sei.cmu.edu/confluence/display/cplusplus/OOP50-CPP.+Do+not+invoke+virtual+functions+from+constructors+or+destructors
   * https://stackoverflow.com/questions/962132/calling-virtual-functions-inside-constructors
   */
  /*
    pre: Assumes global_vertices already initialized
    NOTE have to pass in a lvalue(size_t&) initialzed to 0
  */
  void
  __link_edges_recursive(std::size_t index_x,
                         std::size_t index_y,
                         std::size_t &edge_index,
                         std::vector<bool> &covered);

  /*
    pre: Use only after (edges - vertex) index have built
  */
  void
  __link_faces_recursive(std::size_t index_x,
                         std::size_t index_y,
                         std::size_t &face_index,
                         std::vector<bool> &covered);

public:

  // slice line[a, b] into num of faces which is (num+1) nodes
  /*
    Used only in link_edges_recursive method
  */
  static const std::vector<double>
  __slices_one_dim(double a, double b, std::size_t num);

  /*
    slice num => how many faces in x or y axis:

    form a uniform (m x n) grid with the total of (m*n) faces

    with the out lines as boundary points, and alpha = 0
  */
  simple_square_mesh(double x1 = -1,
                   double x2 = 1,
                   double y1 = -1,
                   double y2 = 1,
                   std::size_t slice_x_num = 10,
                   std::size_t slice_y_num = 10);

  simple_square_mesh(std::pair<double, double> v1,
                   std::pair<double, double> v2,
                   std::size_t slice_x_num = 10,
                   std::size_t slice_y_num = 10);
};

/*
  Mostly copy and paste, which is a very bad style
 */

class simple_triangular_mesh : public simple_mesh {
private:

  std::size_t _x_num;
  std::size_t _y_num;

  void
  __link_edges(size_t v1_index,
               size_t v2_index,
               size_t edge_index);

  void __link_edges_recursive(std::size_t index_x,
                              std::size_t index_y,
                              std::size_t &edge_index,
                              std::vector<bool> &covered);

  void __link_faces_recursive(std::size_t index_x,
                              std::size_t index_y,
                              std::size_t &face_index,
                              std::vector<bool> &covered);

public:
  simple_triangular_mesh(double x1 = -1,
                  double x2 = 1,
                  double y1 = -1,
                  double y2 = 1,
                  std::size_t slice_x_num = 10,
                  std::size_t slice_y_num = 10);

  simple_triangular_mesh(std::pair<double, double> v1,
                  std::pair<double, double> v2,
                  std::size_t slice_x_num = 10,
                  std::size_t slice_y_num = 10);

  const std::vector<face> &
  get_faces() { return _global_faces; };

  const std::vector<vertex> &
  get_vertices() { return _global_vertices; };
};

#endif
