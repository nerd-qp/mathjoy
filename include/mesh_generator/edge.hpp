
/*
  data structure to hold information about an edge
 */

#ifndef _EDGE_
#define _EDGE_

/*
  include stl libraries
 */
#include <memory>
#include <vector>

class vertex;
class face;

//template<std::size_t _vertex_number=2>
class edge {
private:
  std::pair<const vertex *,
            const vertex *> _vertices;

  std::pair<std::size_t,
            std::size_t> _vertices_index;

  std::vector<const face *> _belong_face;
  std::vector<std::size_t> _belong_face_index;
public:
  void add_vertices(const vertex *v1, std::size_t v1_index,
                    const vertex *v2, std::size_t v2_index) {
    _vertices = std::make_pair(v1, v2);
    _vertices_index = std::make_pair(v1_index, v2_index);
  }

  const std::pair<const vertex*,
                  const vertex*> &
  get_vertex() const { return _vertices; };

  const std::pair<std::size_t,
                  std::size_t> &
  get_vertex_index() const { return _vertices_index; };

};

#endif
