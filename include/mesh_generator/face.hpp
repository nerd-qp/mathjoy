
/*
  Data structure to hold information about face
 */

#ifndef _FACE_
#define _FACE_

/*
  Include stl libraries
 */
#include <memory>
#include <vector>
class vertex;
class edge;

class face {
private:
  std::vector<const edge *> _edges;
  std::vector<std::size_t> _edges_index;

public:
  void add_edge(const edge* ed,
                std::size_t ed_index) {
    _edges.emplace_back(ed);
    _edges_index.emplace_back(ed_index);
  }

  const std::vector<const edge *> &
  get_face_edges() { return _edges; };

  std::vector<std::size_t>
  get_face_vertices() const {
    // this is based on the assumption that
    // edges and vertices have the same amount
    // in the same face
    std::vector<std::size_t> res;
    res.reserve(_edges.size());

    for ( const auto &e : _edges ) {
      size_t a, b;

      // use the concept of a priority queue
      std::tie(a, b) = e->get_vertex_index();
      auto lst = {a, b};

      for ( const auto &i : lst ) {
        // if not found, insert
        auto found_iter = std::lower_bound(res.begin(), res.end(), i);
        if ( found_iter == res.end() ||
             *found_iter != i ) {
          res.insert(found_iter, i);
        }
      }
    }

    return res;
  }

};

#endif
