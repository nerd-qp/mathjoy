#include "mesh_generator/mesh_generator.hpp"
#include "utils/utils.hpp"
#include <algorithm>
#include <iostream>

using size_t = std::size_t;

simple_square_mesh::simple_square_mesh(double x1, double x2,
                         double y1, double y2,
                         size_t slice_x_num,
                         size_t slice_y_num)
  : _x_num(slice_x_num+1), _y_num(slice_y_num+1)
{
  if ( slice_x_num <= 0 || slice_y_num <= 0 )
    throw "slices number negative or zero\n";

  _global_vertices = std::vector<vertex>
    (_x_num * _y_num);

  _global_edges = std::vector<edge>
    ((slice_x_num + 1) * slice_y_num + (slice_y_num + 1) * slice_x_num);

  _global_faces = std::vector<face>
    (slice_x_num * slice_y_num);

  auto base_x = __slices_one_dim(x1, x2, slice_x_num);
  auto base_y = __slices_one_dim(y1, y2, slice_y_num);

  std::cout << base_x.size() << ", " << base_y.size() << " slices\n";

  // Construct data structure, link faces, edges, vertices
  for ( std::size_t i = 0; i < base_x.size(); ++ i ) {
    for ( std::size_t j = 0; j < base_y.size(); ++j ) {

      int offset = base_x.size();
      double current_x = base_x[i];
      double current_y = base_y[j];
      vertex *current_vertex = &_global_vertices[i * offset + j];

      // std::cout << current_x << ", " << current_y << " is set\n";

      current_vertex->setX(current_x);
      current_vertex->setY(current_y);

      // Check boundary point
      if ( utils::almost_equal(current_x, x1) ||
           utils::almost_equal(current_x, x2) ||
           utils::almost_equal(current_y, y1) ||
           utils::almost_equal(current_y, y2)) {

        // set default boundary value 0
        current_vertex->setBoundary();

        // std::cout << current_x << ", " << current_y << " is boundary point\n";

      }

    }
  }

  // bool mask to check if the vertex have been traversed
  std::vector<bool> covered(_global_vertices.size(), false);
  size_t init = 0;
  // link edges
  __link_edges_recursive(0, 0, init, covered);
  //  std::cout << init << '\n';

  covered = std::vector<bool>(_global_vertices.size(), false);
  init = 0;
  // link faces
  __link_faces_recursive(0, 0, init, covered);
}

simple_square_mesh::simple_square_mesh(std::pair<double, double> v1,
                         std::pair<double, double> v2,
                         std::size_t slice_x_num,
                         std::size_t slice_y_num) {
  // unpack v1, v2, pass into the last constructor
  double x1, x2, y1, y2;
  std::tie(x1, x2) = v1;
  std::tie(y1, y2) = v2;
  *this = simple_square_mesh(x1, x2, y1, y2, slice_x_num, slice_y_num);
}

void
simple_square_mesh::__link_edges(size_t v1_index,
                          size_t v2_index,
                          size_t edge_index) {
  size_t x = v1_index;
  size_t y = v2_index;
  size_t e = edge_index;

  size_t offset = _x_num;
  vertex *v1 = &_global_vertices[x];
  vertex *v2 = &_global_vertices[y];
  edge *current_edge = &_global_edges[e];

  v1->add_edge(current_edge, e);
  v2->add_edge(current_edge, e);
  current_edge->add_vertices(v1, x, v2, y);

}

void
simple_square_mesh::__link_edges_recursive(size_t index_x,
                                    size_t index_y,
                                    size_t &edge_index,
                                    std::vector<bool> &covered) {
  size_t x = index_x;
  size_t y = index_y;
  size_t offset = _x_num;
  size_t current_index = x * offset + y;

  // if not traversed
  if ( !covered[current_index] ) {

    // if not out of bound, then link
    if ( x + 1 < _x_num ) {
      size_t right_index = (x+1) * offset + y;

      // increment edge index here
      __link_edges(right_index, current_index, edge_index++);
      __link_edges_recursive(x+1, y, edge_index, covered);
    }

    // if not out of bound, then link
    if ( y + 1 < _y_num ) {
      size_t upper_index = x * offset + y + 1;

      // increment edge index here
      __link_edges(upper_index, current_index, edge_index++);
      __link_edges_recursive(x, y+1, edge_index, covered);
    }
    covered[current_index] = true;
  }
}

void
simple_square_mesh::__link_faces_recursive(size_t index_x,
                                    size_t index_y,
                                    size_t &face_index,
                                    std::vector<bool> &covered) {
  size_t x = index_x;
  size_t y = index_y;
  size_t offset = _x_num;
  size_t current_index = x * offset + y;

  //std::cout << '(' << x << ',' << y << ')' << '\n';
  //std::cout << face_index << '\n';

  // check if out of bound
  // on the diagonal side

  if ( !covered[current_index] ) {
    if ( x+1 < _x_num && y+1 < _y_num ) {

      using con_iter = std::vector<vertex>::const_iterator;

      // link edges to face
      // find vectices of the region
      con_iter v1 = _global_vertices.cbegin() + current_index;
      con_iter v2 = _global_vertices.cbegin() + current_index + offset;
      con_iter v3 = _global_vertices.cbegin() + current_index + offset + 1;
      con_iter v4 = _global_vertices.cbegin() + current_index + 1;

      // form a cyclic tuple for std::accumulate
      auto local_vs = {v1, v2, v3, v4};
      auto local_vs_shift = {v2, v3, v4, v1};

      // find common edges
      std::vector<size_t> common_edges;
      utils::zip(local_vs.begin(), local_vs.end(),
                 local_vs_shift.begin(), local_vs_shift.end(),
                 [&] (con_iter a, con_iter b) {
                   auto a_edges = a->getEdges();
                   auto b_edges = b->getEdges();

                   size_t result =
                     *std::find_first_of(a_edges.begin(), a_edges.end(),
                                         b_edges.begin(), b_edges.end());
                   common_edges.emplace_back(result);
                 });

      for ( const auto &e : common_edges )
        _global_faces[face_index].add_edge(&_global_edges[e], e);
      //   std::cout << e << ' ';
      // std::cout << '\n';


      // increment face index after linking
      face_index++;

      // recursive calls
      __link_faces_recursive(x+1, y, face_index, covered);
      __link_faces_recursive(x, y+1, face_index, covered);

    }

    covered[current_index] = true;
  }

}


const std::vector<double>
simple_square_mesh::__slices_one_dim(double a, double b, std::size_t num) {
  //
  if ( num == 1 )
    return {a, b};

  double length = (b - a) / static_cast<double>(num);

  std::vector<double> res_vertex(num+1);

  for ( std::size_t i = 0; i < res_vertex.size(); ++ i ) {
    res_vertex[i] = a + i * length;
  }

  return res_vertex;
}

void
simple_triangular_mesh::__link_edges(size_t v1_index,
                                     size_t v2_index,
                                     size_t edge_index) {
  size_t x = v1_index;
  size_t y = v2_index;
  size_t e = edge_index;

  size_t offset = _x_num;
  vertex *v1 = &_global_vertices[x];
  vertex *v2 = &_global_vertices[y];
  edge *current_edge = &_global_edges[e];

  v1->add_edge(current_edge, e);
  v2->add_edge(current_edge, e);
  current_edge->add_vertices(v1, x, v2, y);

}


void
simple_triangular_mesh::__link_edges_recursive(size_t index_x,
                                               size_t index_y,
                                               size_t &edge_index,
                                               std::vector<bool> &covered) {
  size_t x = index_x;
  size_t y = index_y;
  size_t offset = _x_num;
  size_t current_index = x * offset + y;

  // if not traversed
  if ( !covered[current_index] ) {

    // link upper right vertex
    if ( x + 1 < _x_num && y + 1 < _y_num ) {
      size_t upper_right_index = (x+1) * offset + y + 1;

      __link_edges(upper_right_index, current_index, edge_index++);
      // no need to traverse
    }

    // if not out of bound, then link
    if ( x + 1 < _x_num ) {
      size_t right_index = (x+1) * offset + y;

      // increment edge index here
      __link_edges(right_index, current_index, edge_index++);
      __link_edges_recursive(x+1, y, edge_index, covered);
    }

    // if not out of bound, then link
    if ( y + 1 < _y_num ) {
      size_t upper_index = x * offset + y + 1;

      // increment edge index here
      __link_edges(upper_index, current_index, edge_index++);
      __link_edges_recursive(x, y+1, edge_index, covered);
    }

    covered[current_index] = true;
  }
}

void simple_triangular_mesh::
__link_faces_recursive(size_t index_x,
                       size_t index_y,
                       size_t &face_index,
                       std::vector<bool> &covered) {
  size_t x = index_x;
  size_t y = index_y;
  size_t offset = _x_num;
  size_t current_index = x * offset + y;

  //std::cout << '(' << x << ',' << y << ')' << '\n';
  //std::cout << face_index << '\n';

  // check if out of bound
  // on the diagonal side

  if ( !covered[current_index] ) {
    if ( x+1 < _x_num && y+1 < _y_num ) {

      using con_iter = std::vector<vertex>::const_iterator;

      // link edges to face
      // find vectices of the region
      con_iter v1 = _global_vertices.cbegin() + current_index;
      con_iter v2 = _global_vertices.cbegin() + current_index + offset;
      con_iter v3 = _global_vertices.cbegin() + current_index + offset + 1;
      con_iter v4 = _global_vertices.cbegin() + current_index + 1;

      // form a cyclic tuple for std::accumulate
      std::vector<con_iter> local_vs = {v1, v2, v3};
      std::vector<con_iter> local_vs_shift = {v2, v3, v1};

      // find common edges
      std::vector<size_t> common_edges;
      utils::zip(local_vs.begin(), local_vs.end(),
                 local_vs_shift.begin(), local_vs_shift.end(),
                 [&] (con_iter a, con_iter b) {
                   auto a_edges = a->getEdges();
                   auto b_edges = b->getEdges();

                   size_t result =
                     *std::find_first_of(a_edges.begin(), a_edges.end(),
                                         b_edges.begin(), b_edges.end());
                   common_edges.emplace_back(result);
                 });

      for ( const auto &e : common_edges ) {
        //std::cout << "edge index: " << e << ' ';
        //std::cout << e << ' ';
        _global_faces[face_index].add_edge(&_global_edges[e], e);
      }
      //std::cout << '\n';

      // increment index
      ++ face_index;

      // clear it
      common_edges.clear();

      // don't do this mate...
      // NOTE std::initialzer_list is const
      local_vs = {v1, v4, v3};
      local_vs_shift = {v4, v3, v1};

      // fill again
      utils::zip(local_vs.begin(), local_vs.end(),
                 local_vs_shift.begin(), local_vs_shift.end(),
                 [&] (con_iter a, con_iter b) {
                   auto a_edges = a->getEdges();
                   auto b_edges = b->getEdges();

                   size_t result =
                     *std::find_first_of(a_edges.begin(), a_edges.end(),
                                         b_edges.begin(), b_edges.end());
                   common_edges.emplace_back(result);
                 });


      for ( const auto &e : common_edges ) {
        //        std::cout << e << ' ';
        _global_faces[face_index].add_edge(&_global_edges[e], e);
      }
      //std::cout << '\n';

      // increment face index after linking
      ++ face_index;

      // recursive calls
      __link_faces_recursive(x+1, y, face_index, covered);
      __link_faces_recursive(x, y+1, face_index, covered);

    }

    covered[current_index] = true;
  }

}

simple_triangular_mesh::
simple_triangular_mesh(double x1, double x2,
                       double y1, double y2,
                       size_t slice_x_num,
                       size_t slice_y_num)
  : _x_num(slice_x_num+1), _y_num(slice_y_num+1)
{
  if ( slice_x_num <= 0 || slice_y_num <= 0 )
    throw "slices number negative or zero\n";

  _global_vertices = std::vector<vertex>
    (_x_num * _y_num);

  _global_edges = std::vector<edge>
    (3 * slice_y_num * slice_x_num + slice_x_num + slice_y_num);
    // ((slice_x_num + 1) * slice_y_num +
    //  (slice_y_num + 1) * slice_x_num +
    //  slice_x_num * slice_y_num);

  _global_faces = std::vector<face>
    (slice_x_num * slice_y_num * 2);

  auto base_x = simple_square_mesh::__slices_one_dim(x1, x2, slice_x_num);
  auto base_y = simple_square_mesh::__slices_one_dim(y1, y2, slice_y_num);

  std::cout << base_x.size() << ", " << base_y.size() << " slices\n";

  // Construct data structure, link faces, edges, vertices
  for ( std::size_t i = 0; i < base_x.size(); ++ i ) {
    for ( std::size_t j = 0; j < base_y.size(); ++j ) {

      int offset = base_x.size();
      double current_x = base_x[i];
      double current_y = base_y[j];
      vertex *current_vertex = &_global_vertices[i * offset + j];

      // std::cout << current_x << ", " << current_y << " is set\n";

      current_vertex->setX(current_x);
      current_vertex->setY(current_y);

      // Check boundary point
      if ( utils::almost_equal(current_x, x1) ||
           utils::almost_equal(current_x, x2) ||
           utils::almost_equal(current_y, y1) ||
           utils::almost_equal(current_y, y2)) {

        // set default boundary value 0
        current_vertex->setBoundary();

        // std::cout << current_x << ", " << current_y << " is boundary point\n";

      }

    }
  }

  // bool mask to check if the vertex have been traversed
  std::vector<bool> covered(_global_vertices.size(), false);
  size_t init = 0;
  // link edges
  __link_edges_recursive(0, 0, init, covered);
  //  std::cout << init << '\n';

  covered = std::vector<bool>(_global_vertices.size(), false);
  init = 0;
  // link faces
  __link_faces_recursive(0, 0, init, covered);
}

simple_triangular_mesh::
simple_triangular_mesh(std::pair<double, double> v1,
                       std::pair<double, double> v2,
                       size_t slice_x_num,
                       size_t slice_y_num)
{ // unpack v1, v2, pass into the last constructor
  double x1, x2, y1, y2;
  std::tie(x1, x2) = v1;
  std::tie(y1, y2) = v2;
  *this = simple_triangular_mesh(x1, x2, y1, y2, slice_x_num, slice_y_num);
}
