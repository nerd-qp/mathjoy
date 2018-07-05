// included due to numeric.h didn't include it

#include "FVM_2dim.hpp"

// for benchmarking
#include <chrono>
// for make_pair
#include <utility>

namespace FVM {
  void FVM_solver_2dim::
  _connect(std::vector<std::vector<std::size_t>>& triangles_indices,
           std::vector<std::vector<std::size_t>>& self_indices,
           std::vector<std::vector
           <std::pair<std::size_t, std::size_t>>>& other_indices) {
      auto size = _points.size();
      // maybe unecessary
      triangles_indices.clear();
      self_indices.clear();
      other_indices.clear();
      triangles_indices.resize(size);
      self_indices.resize(size);
      other_indices.resize(size);

      std::size_t index = 0;
      // assume the index given will not be out of bound
      for ( const auto& triangle : _triangulation.getTriangles() ) {
        const auto& index1 = triangle.getIndex1();
        const auto& index2 = triangle.getIndex2();
        const auto& index3 = triangle.getIndex3();
        triangles_indices[index1].emplace_back(index);
        triangles_indices[index2].emplace_back(index);
        triangles_indices[index3].emplace_back(index);
        self_indices[index1].emplace_back(1);
        self_indices[index2].emplace_back(2);
        self_indices[index3].emplace_back(3);

        other_indices[index1].emplace_back(std::make_pair(index2, index3));
        other_indices[index2].emplace_back(std::make_pair(index1, index3));
        other_indices[index3].emplace_back(std::make_pair(index1, index2));
        index++;
      }

  }


  // void FVM_solver_2dim::
  // _connect(std::vector<std::vector<
  //          std::pair<std::size_t, std::size_t>>>& connect_mapping) {
  //   auto size = _points.size();
  //   // maybe unecessary
  //   connect_mapping.clear();
  //   connect_mapping.resize(size);

  //   std::size_t index = 0;
  //   // assume the index given will not be out of bound
  //   for ( const auto& triangle : _triangulation.getTriangles() ) {
  //     const auto& index1 = triangle.getIndex1();
  //     const auto& index2 = triangle.getIndex2();
  //     const auto& index3 = triangle.getIndex3();
  //     connect_mapping[index1].emplace_back(std::make_pair(index, 1));
  //     connect_mapping[index2].emplace_back(std::make_pair(index, 2));
  //     connect_mapping[index3].emplace_back(std::make_pair(index, 3));
  //     index++;
  //   }

  //   // for ( std::size_t i = 0u; i < size; ++ i ) {

  //     // // wrong algo, could possibly have some common edge, adding redundant vertex
  //     // // find which triangle contains this point
  //     // // use the almost equal method provided by delaunary triangulation
  //     // // and add the other 2 points to the connect_mapping
  //     // auto cbegin = _triangulation.getTriangles().cbegin();
  //     // auto cend   = _triangulation.getTriangles().cend();
  //     // std::size_t index = 0;

  //     // std::for_each(cbegin, cend,
  //     //               [&connect_mapping, i, this, &index]
  //     //               (const Triangle<double>& triangle) {
  //     //                 const auto& point = _points[i];

  //     //                 // cast to int
  //     //                 int mask[3] = {almost_equal(point, triangle.p1),
  //     //                                almost_equal(point, triangle.p2),
  //     //                                almost_equal(point, triangle.p3)};

  //     //                 int choice = mask[0] + (mask[1] << 1) + (mask[2] << 2);
  //     //                 // unspecified
  //     //                 auto _the_pair = std::make_pair(index, 0ul);
  //     //                 // see which point it fits, if none, then move on
  //     //                 switch( choice ) {
  //     //                   // the possibly the common case
  //     //                 case 0: {
  //     //                   break;
  //     //                 }
  //     //                 case 1: {
  //     //                   _the_pair.second = 1;
  //     //                   connect_mapping[i].push_back(_the_pair);
  //     //                   connect_mapping[i].push_back(_the_pair);
  //     //                   break;
  //     //                 }
  //     //                 case 2: {
  //     //                   _the_pair.second = 2;
  //     //                   connect_mapping[i].push_back(_the_pair);
  //     //                   connect_mapping[i].push_back(_the_pair);
  //     //                   break;
  //     //                 }
  //     //                 case 4: {
  //     //                   _the_pair.second = 3;
  //     //                   connect_mapping[i].push_back(_the_pair);
  //     //                   connect_mapping[i].push_back(_the_pair);
  //     //                   break;
  //     //                 }
  //     //                 default: {
  //     //                   // no one should get here, just throw
  //     //                   throw "some thing wrong with the triangle";
  //     //                 }
  //     //                 }
  //     //                 // increment index
  //     //                 ++index;
  //     //               });
  //     // // correct algo
  //     // // find which edge contains this point
  //     // // and add the other vertex into the list
  //     // // because every edge has two vertices, no redundant neighbour could be added
  //     // // NOTE: this code is highly inefficient
  //     // auto cbegin = _triangulation.getEdges().cbegin();
  //     // auto cend   = _triangulation.getEdges().cend();

  //     // std::for_each(cbegin, cend,
  //     //               [&connect_mapping, i, this]
  //     //               (const Edge<double>& edge) {
  //     //                 const auto& point = _points[i];

  //     //                 // cast to int
  //     //                 int mask[2] = {almost_equal(_points[i], edge.p1),
  //     //                                almost_equal(_points[i], edge.p2)};

  //     //                 int choice = mask[0] + (mask[1] << 1);

  //     //                 // see which point it fits, if none, then move on
  //     //                 switch( choice ) {
  //     //                   // the possibly the common case
  //     //                 case 0: {
  //     //                   break;
  //     //                 }
  //     //                 case 1: {
  //     //                   connect_mapping[i].push_back(edge.p2);
  //     //                   break;
  //     //                 }
  //     //                 case 2: {
  //     //                   connect_mapping[i].push_back(edge.p1);
  //     //                   break;
  //     //                 }
  //     //                 default: {
  //     //                   // no one should get here, just throw
  //     //                   throw "some thing wrong with the triangle";
  //     //                 }
  //     //                 }
  //     //               });
  //   // }
  // }

  void FVM_solver_2dim::_circumCenter
  (// the triangle indices should not be const
   // so it can be reArrange into neighbour
   std::vector<std::vector<std::size_t>>& triangles_indices,
   std::vector<std::vector<std::size_t>>& self_indices,
   std::vector<std::vector
   <std::pair<std::size_t, std::size_t>>>& other_indices,
   std::vector<std::vector<Vector2<double>>>& circumCenter_points) {

    // clear and reserve space
    circumCenter_points.clear();
    auto size = triangles_indices.size();
    circumCenter_points.resize(size);
    for ( std::size_t i = 0; i < size; ++ i )
      circumCenter_points[i].reserve(triangles_indices[i].size());

    // Here assume that every point that's not on the boundary will have a fully
    // connected edges meaning that (n triangles == n neighbour) sometimes would
    // ignore the actual shape of the triangulation

    // compute circumcenter of each triangle(maybe the original shape) and put
    // into data structure
    const auto& _triangles = _triangulation.getTriangles();
    for ( std::size_t i = 0; i < size; ++ i ) {
      // triangle indices of node i, non const for arrange
      auto& local_triangles_indices = triangles_indices[i];
      // const for read only
      auto& local_other_indices = other_indices[i];
      auto& local_self_indices = self_indices[i];

      // if the point is on boundary, no need to compute(or arrange)
      if ( _on_boundary[i] ) { continue; }

      // triangles.reserve(size);
      // self.reserve(size);
      // for ( const auto& connect : connect_mapping ) {
      //   triangles.emplace_back(connect.first);
      //   self.emplace_back(connect.second);
      // }

      // extract information on two other vertices
      // std::vector<std::pair<std::size_t, std::size_t>> other_indices;
      //_markOtherIndices(self_index, other_indices);

      // reArrange the indices so that the triangles with the common edge is side
      // by side
      _reArrange(local_triangles_indices, local_self_indices,
                 local_other_indices);

      // now we can compute the circumcenter
      for ( std::size_t j = 0; j < local_triangles_indices.size(); ++ j ) {
        const auto& local_triangle = _triangles[local_triangles_indices[j]];
        circumCenter_points[i].emplace_back(local_triangle.circumCenter());
      }
    }
  }

  // void FVM_solver_2dim::
  // _markOtherIndices(const std::vector<std::size_t>& self_index,
  //                   <std::pair<std::size_t, std::size_t>>& other_indices) {

  //   const auto& triangles = _triangulation.getTriangles();
  //   const auto& size = triangles.size();
  //   other_indices.clear();
  //   other_indices.resize(size);
  //   for ( std::size_t i = 0; i < size; ++ i ) {
  //     switch(self_index[i]) {
  //     case 1: {
  //       other_indices[i] = std::make_pair(triangles[i].index2, triangles[i].index3);
  //       break;
  //     }
  //     case 2: {
  //       other_indices[i] = std::make_pair(triangles[i].index1, triangles[i].index3);
  //       break;
  //     }
  //     case 3: {
  //       other_indices[i] = std::make_pair(triangles[i].index1, triangles[i].index2);
  //       break;
  //     }
  //     default: {
  //       throw "no one should get here!\n";
  //     }
  //     }
  //   }

  //   return;
  // }

  void FVM_solver_2dim::
  _reArrange(std::vector<std::size_t>& triangles,
             std::vector<std::size_t>& self_indices,
             std::vector
             <std::pair<std::size_t, std::size_t>>& other_indices) {

    auto size = triangles.size();

    // if the node is linked to only 2 triangles, not possible to have link head
    // and tail, which means useless to perform reArrange
    if ( size <= 2 ) {
      std::cerr << "we have some not so well connected triangles"
                << "(only 2 or less triangle)\n";
      return;
    }

    // run find iteratively, and swap(rotate to neighbour)
    for ( std::size_t i = 0; i < size-1; ++ i ) {
      // keep j here for checking afterwards
      std::size_t j = i+1;
      // find the first neighbour
      for ( ; j < size; ++ j ) {
        // the four other indices(two for each triangle)
        const auto& first1 = other_indices[i].first;
        const auto& first2 = other_indices[j].first;
        const auto& second1 = other_indices[i].second;
        const auto& second2 = other_indices[j].second;

        bool mask[4] = {first1 == first2,
                        first1 == second2,
                        second1 == first2,
                        second1 == second2};

        auto sum = mask[0] + (mask[1] << 1) + (mask[2] << 2) + (mask[3] << 3);

        // here assume it's not possible that two triangles are the same
        // therefore the four cases would only have one possible case each time
        // here(out of lazyness) I treat self_indices as the P
        switch ( sum ) {
        case 0: {
          break;
        }
        case 1: {
          self_indices[i] = first1;
          break;
        }
        case 2: {
          self_indices[i] = first1;
          break;
        }
        case 4: {
          self_indices[i] = second1;
          break;
        }
        case 8: {
          self_indices[i] = second1;
          break;
        }
        default: {
          std::cerr << first1 << ' '
                    << first2 << ' '
                    << second1 << ' '
                    << second2 << '\n';
          std::cerr << i << ' ' << j << '\n';
          throw "shouldn't be here though\n";
        }
        }

        // if have a common edge, neighbour found
        if ( sum != 0 ) {
          // swap it back
          if ( i+1 != j ) {
            std::swap(triangles[i + 1], triangles[j]);
            std::swap(other_indices[i + 1], other_indices[j]);
            //std::swap(self_indices[i + 1], self_indices[j]);
          }

          // stop the search, found the first one
          // continue on the next iteration
          break;
        }
      }

      // if we don't find one... that means the triangulation is bad...
      if ( j == size ) {
        // std::cerr << first1 << ' '
        //           << first2 << ' '
        //           << second1 << ' '
        //           << second2 << '\n';
        // bad thingy
        std::cerr << i << ' ';
        std::cerr << j << ' ';
        std::cerr << "we have some not so well connected triangles\n";
        return;
      }
    }
    // try to link the head and tail if the common node isn't already linked
    // and if it's already linked, that means it's a bad triangulation with only
    // two triangles around this node
    // in such case, just check the size of the triangulation, this case is
    // case out by the if above
    auto head = 0;
    auto tail = other_indices.size() - 1;
    const auto& first1 = other_indices[head].first;
    const auto& first2 = other_indices[tail].first;
    const auto& second1 = other_indices[head].second;
    const auto& second2 = other_indices[tail].second;
    bool mask[4] = {first1 == first2,
                    first1 == second2,
                    second1 == first2,
                    second1 == second2};

    auto sum = mask[0] + (mask[1] << 1) + (mask[2] << 2) + (mask[3] << 3);

    switch ( sum ) {
    case 0: {
      // cannot link the head and the tail
      std::cerr << "we have some not so well connected triangles"
                << "(cannot connect head and tail)\n";
      return;
      //      break;
    }
    case 1: {
      self_indices[tail] = first1;
    }
    case 2: {
      self_indices[tail] = first1;
      break;
    }
    case 4: {
      self_indices[tail] = second1;
      break;
    }
    case 8: {
      self_indices[tail] = second1;
      break;
    }
    default: {
      std::cerr << first1 << ' '
                << first2 << ' '
                << second1 << ' '
                << second2 << '\n';
      throw "shouldn't be here though\n";
    }
    }
  }

  std::vector<Edge<double>> FVM_solver_2dim::
  assemble_all() {
    // connectivity information of the points
    // index i connects to a vector of index
    // std::vector<std::vector<Vector2<double>>> connect_mapping;
    // std::vector<std::vector<std::pair<std::size_t, std::size_t>>> connect_mapping;

    // information on triangle index and position of the node in the triangle
    std::vector<std::vector<std::size_t>> triangles_indices;
    std::vector<std::vector<std::size_t>> self_indices;
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> other_indices;

    // each correspong to one triangle
    std::vector<std::vector<Vector2<double>>> circumCenter_points;

    // get the connectivity information
    auto start = std::chrono::high_resolution_clock::now();

    this->_connect(triangles_indices, self_indices, other_indices);
    this->_circumCenter(triangles_indices, self_indices,
                        other_indices, circumCenter_points);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "takes " << diff.count() << " to complete connecting\n";

    auto size = other_indices.size();

    // showing the components
    for ( auto i = 0u; i < size; ++ i ) {

      // don't show point on boundary
      if ( _on_boundary[i] ) {
        continue;
      }

      std::cout << _points[i] << " => ";
      auto index = 0u;
      for ( const auto& neighbour : other_indices[i] ) {
        std::cout << '<' << neighbour.first << ", " << neighbour.second << '>' << ' ';
        // the Pi+1
        std::cout << self_indices[i][index] << ' ';
        // std::cout << circumCenter_points[i][index] << ' ';
        ++index;
        // std::cout << neighbour << ' ';
      }

      std::cout << size << '\n';
    }

    for ( const auto& edge : _triangulation.getEdges() )  {
      //std::cout << edge.index1 << ' ' << edge.index2 << '\n';
    }

    // now we can assemble the stiffnees matrix
    std::vector<Edge<double>> edges;

    for ( std::size_t i = 0u; i < size; ++ i ) {

      if ( _on_boundary[i] ) {
        // if the node is on boundary
        // simply assign...
        // so the weight should be 1
        // the right side should be alpha
        _stiffness_matrix(i, i) = 1;
        auto alpha = 0.0;
        _right_f_phi(i) = alpha;
        continue;
      }

      const auto& points = circumCenter_points[i];
      auto size2 = points.size();

      // sum of all weights for U(P_0)
      double sum = 0;
      // the right side integral
      double right_sum = 0;
      const auto& PB = _points[i];

      // assemble the stiffness matrix
      for ( std::size_t j = 0u; j < size2; ++ j ) {
        // q_i, q_i+1
        const auto& cirPointA = points[j];
        const auto& cirPointB = points[(j+1) % size2];

        // P_0, P_i+1
        const auto& PA = _points[self_indices[i][j]];

        auto numerator = cirPointA.dist(cirPointB);
        auto denominator = PB.dist(PA);

        auto weight = numerator / denominator;
        // std::cout << weight << '\n';

        // assemble on weights to P_i
        // _stiffness_matrix(i, self_indices[i][j]) += -weight;
        _stiffness_matrix(i, self_indices[i][j]) = -weight;

        // accumulate
        sum += weight;

        right_sum +=
          _quadrature_table->quadrature_compute
          (two_dimension::canonical::
           function2Canonical(_f, std::make_pair(PB.x, PB.y),
                              std::make_pair(cirPointA.x, cirPointA.y),
                              std::make_pair(cirPointB.x, cirPointB.y)));

        // create edge connection, head of point would connect to tail of point
        edges.emplace_back(cirPointA, cirPointB, -1, -1);
        // edges.emplace_back(cirPointA, PB, -1, -1);
      }
      // assemble to P_0
      // _stiffness_matrix(i, i) += sum;
      _stiffness_matrix(i, i) = sum;
      _right_f_phi(i) = right_sum;
    }

    return edges;
  }
}
