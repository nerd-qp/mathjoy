
/*
  Holding information about 2 dimensional vertices

  Store additional information about containing faces
 */

#ifndef _VERTEX_
#define _VERTEX_

/*
  Include stl libraries
 */
#include <memory>
#include <vector>

// forward declaration
class edge;
class face;

// vertex class
class vertex {
private:
  double                              _x;
  double                              _y;
  std::vector<const edge *>           _belong_edge;
  std::vector<std::size_t>            _belong_edge_index;
  std::pair<bool, double>             _boundary_info;

public:
  vertex() : _boundary_info(std::make_pair(false, 0)) {};
  void setX(double x) { _x = x; };
  void setY(double y) { _y = y; };
  void setBoundary(double alpha = 0) {
    _boundary_info = std::make_pair(true, alpha);
  }

  double getX() const { return _x; };
  double getY() const { return _y; };
  std::pair<bool, double> getBoundary() const {
    return _boundary_info;
  }

  void add_edge(const edge * ref_edge, std::size_t ref_index) {
    _belong_edge.emplace_back(ref_edge);
    _belong_edge_index.emplace_back(ref_index);
  }

  const std::vector<std::size_t>&
  getEdges() const { return _belong_edge_index; };

};

#endif
