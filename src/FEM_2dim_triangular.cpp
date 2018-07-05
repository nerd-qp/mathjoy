
/*
 * Implementation of 2 dimension FEM solver
 */

/*
 * Include libraries
 */
#include "FEM_2dim.hpp"
#include <chrono>
#include <cmath>

// for testing only
#include <boost/numeric/ublas/io.hpp>

// progress display
#include <boost/progress.hpp>

namespace linear {

  /*
   * A linear coordinate transformation
   * Take two vector of x, y coordinates
   * return a pair of function that maps kexi, yita
   * from [-1, 1] x [-1, 1] to the corresponding region
   * represented by x, y coordinates
   */
  inline std::pair<utils::function, utils::function>
  FEM_triangular::__linear_co_transform(const vector<double> &x_coor,
                                        const vector<double> &y_coor) const {
    if ( x_coor.size() != 3 || x_coor.size() != y_coor.size() ) {
      throw "coordinates size not qualified!\n";
    }

    double x_i = x_coor[0], x_j = x_coor[1], x_k = x_coor[2];
    double y_i = y_coor[0], y_j = y_coor[1], y_k = y_coor[2];
    utils::function To_x([=] (double kexi, double yita) -> double {
                          return x_i + (x_j - x_i) * kexi + (x_k - x_i) * yita;
                        });

    utils::function To_y([=] (double kexi, double yita) -> double {
                          return y_i + (y_j - y_i) * kexi + (y_k - y_i) * yita;
                        });

    return std::make_pair(std::move(To_x), std::move(To_y));
  }

  /*
   * Compute the derivative function for x, y coordinates
   */
  inline std::tuple<utils::function,
                    utils::function,
                    utils::function,
                    utils::function>
  FEM_triangular::__local_derivative_function_get
  (const vector<double> &x_coor,
   const vector<double> &y_coor) const {
    if ( x_coor.size() != 3 || x_coor.size() != y_coor.size() ) {
      throw "coordinates size not qualified!\n";
    }

    utils::function x_kexi, x_yita, y_kexi, y_yita;

    double x_i = x_coor[0], x_j = x_coor[1], x_k = x_coor[2];
    double y_i = y_coor[0], y_j = y_coor[1], y_k = y_coor[2];

    x_kexi = ([=] (double kexi, double yita) -> double {
                return x_j - x_i;
              });

    x_yita = ([=] (double kexi, double yita) -> double {
                return x_k - x_i;
              });

    y_kexi = ([=] (double kexi, double yita) -> double {
                return y_j - y_i;
              });

    y_yita = ([=] (double kexi, double yita) -> double {
                return y_k - y_i;
              });

    return std::make_tuple(std::move(x_kexi),
                           std::move(x_yita),
                           std::move(y_kexi),
                           std::move(y_yita));
  }

  /*
   * Some helper funciton to get the corresponding
   * x_coor, y_coor and vertices_index
   * since multiple places involved
   * Please forgive me for the name...
   * personal use only...
   */

  inline std::tuple<std::vector<size_t>, vector<double>, vector<double>>
  FEM_triangular::__index_x_coor_y_coor_get(const face &element_face) const {
    auto vertices_index = element_face.get_face_vertices();
    size_t size = vertices_index.size();
    vector<double> x_coor(size);
    vector<double> y_coor(size);

    const std::vector<vertex> &global_ver = _mesh->get_vertices();

    for ( size_t i = 0; i < size; ++ i) {
      size_t index = vertices_index[i];
      x_coor[i] = global_ver[index].getX();
      y_coor[i] = global_ver[index].getY();
    }

    return std::make_tuple(vertices_index, x_coor, y_coor);
  }

  /*
   * As suggested by the name
   * Pre: Have to prepare the mesh and have generated all
   * required base functions
   */
  inline void
  FEM_triangular::__assemble_local_stiffness_matrix(const face &element_face) {
    // simply copy and paste here

    // prepare data
    matrix<utils::function> K_fun_matrix, M_fun_matrix;
    matrix<std::pair<size_t, size_t>> K_global_index_matrix, M_global_index_matrix;
    matrix<double> K_matrix, M_matrix;

    // auto pre = std::chrono::high_resolution_clock::now();

    // generate the local functions
    std::tie(K_fun_matrix, K_global_index_matrix) =
      K_element_matrix(element_face);

    std::tie(M_fun_matrix, M_global_index_matrix) =
      M_element_matrix(element_face);
    // std::cout << "Time used to initialize: " <<
    //   std::chrono::duration_cast<std::chrono::duration<double>>
    //   (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;

    // pre = std::chrono::high_resolution_clock::now();

    // compute the local matrix
    K_matrix = _quadrature_table->quadrature_compute(K_fun_matrix);
    M_matrix = _quadrature_table->quadrature_compute(M_fun_matrix);

    // std::cout << "Time used to compute: " <<
    //   std::chrono::duration_cast<std::chrono::duration<double>>
    //   (std::chrono::high_resolution_clock::now() - pre).count() << std::endl;

    size_t size1 = K_fun_matrix.size1();
    size_t size2 = K_fun_matrix.size2();
    size_t first_index, second_index;

    // assemble global matrix
    for ( size_t i = 0; i < size1; ++ i ) {
      for ( size_t j = 0; j < size2; ++ j ) {
        auto k_index = &K_global_index_matrix(i, j);
        auto m_index = &M_global_index_matrix(i, j);

        first_index = k_index->first;
        second_index = k_index->second;
        _stiffness_matrix (first_index, second_index) += K_matrix(i, j);

        first_index = m_index->first;
        second_index = m_index->second;

        _stiffness_matrix (first_index, second_index) += M_matrix(i, j);
      }
    }


  }

  inline void
  FEM_triangular::__assemble_local_right_f_phi(const face &element_face) {
    vector<utils::function> f_fun_vector;
    vector<size_t> f_global_index_vector;
    vector<double> f_matrix;

    std::tie(f_fun_vector, f_global_index_vector) =
      f_element_vector(element_face);

    f_matrix = _quadrature_table->quadrature_compute(f_fun_vector);

    size_t size = f_fun_vector.size();
    size_t index;

    for ( size_t i = 0; i < size; ++ i ) {
      auto f_index = &f_global_index_vector(i);

      index = *f_index;
      _right_f_phi (index) += f_matrix( i );
    }

  }

  /*
   * delegating constructor
   */

  FEM_triangular::FEM_triangular(std::pair<double, double> v1,
                                 std::pair<double, double> v2,
                                 size_t slice_x_num,
                                 size_t slice_y_num,
                                 utils::function p,
                                 utils::function q,
                                 utils::function f,
                                 size_t degree,
                                 const size_t &freedom) :
    FEM_solver_2dim_base(p, q, f, (slice_x_num + 1) * (slice_y_num + 1)),
    _base_fun_vec(freedom),
    _base_diff_vec_kexi(freedom),
    _base_diff_vec_yita(freedom) {
    // should be fixed and not upgradable,
    // since the generation of shape function is hard to write
    if ( degree > 1 )
      throw "only support degree of 1 now!\n";
    _mesh = utils::make_unique<simple_triangular_mesh>(v1, v2,
                                                       slice_x_num,
                                                       slice_y_num);

    _quadrature_table =
      utils::make_unique<two_dimension::canonical::triangular>(5);

    _base_fun_vec.resize(3);

    _base_fun_vec[0] = utils::function([] (double x, double y) -> double {
                                         return 1 - x - y;
                                       });

    _base_fun_vec[1] = utils::function([] (double x, double y) -> double {
                                         return x;
                                       });

    _base_fun_vec[2] = utils::function([] (double x, double y) -> double {
                                         return y;
                                       });

    _base_diff_vec_kexi.resize(3);

    _base_diff_vec_kexi[0] = utils::function([] (double x, double y) -> double {
                                               return -1;
                                             });

    _base_diff_vec_kexi[1] = utils::function([] (double x, double y) -> double {
                                               return 1;
                                             });

    _base_diff_vec_kexi[2] = utils::function([] (double x, double y) -> double {
                                               return 0;
                                             });

    _base_diff_vec_yita.resize(3);

    _base_diff_vec_yita[0] = utils::function([] (double x, double y) -> double {
                                               return -1;
                                             });

    _base_diff_vec_yita[1] = utils::function([] (double x, double y) -> double {
                                               return 0;
                                             });

    _base_diff_vec_yita[2] = utils::function([] (double x, double y) -> double {
                                               return 1;
                                             });

    assemble_stiffness_matrix();
    assemble_right_f_phi();

    for ( size_t i = 0; i < _mesh->get_vertices().size(); ++ i )
      handle_boundary_point(i);

  }

  /*
   * Pre: mesh generation complete
   */
  void
  FEM_triangular::assemble_stiffness_matrix() {
    auto all_faces = _mesh->get_faces();

    // for ( const auto &face : all_faces ) {
    //   __assemble_local_stiffness_matrix(face);
    //   std::cout << '.';
    // }

    boost::progress_display show_progress(all_faces.size());

    for ( size_t i = 0; i < all_faces.size(); ++ i ) {

      __assemble_local_stiffness_matrix(all_faces[i]);

      ++show_progress;
      // Given the optimization flag, have to use flush to flush output
      // if ( (i * 100 / all_faces.size()) % 10 == 0 )
      //   std::cout << '.' << std::flush;
    }

  }

  void
  FEM_triangular::assemble_right_f_phi() {
    auto all_faces = _mesh->get_faces();

    boost::progress_display show_progress(all_faces.size());
    for ( size_t i = 0; i < all_faces.size(); ++ i ) {
      __assemble_local_right_f_phi(all_faces[i]);

      // have to use flush
      // if ( (i * 100 / all_faces.size()) % 10 == 0 )
      //   std::cout << i << '.' << std::flush;

      ++show_progress;
    }
  }

  /*
   * Should only be used post stiffness matrix assembling
   * and right_f_phi assembing done
   */
  inline void
  FEM_triangular::handle_boundary_point(size_t index) {
    bool isBoundary;
    double alpha;
    auto v = &(_mesh->get_vertices()[index]);
    std::tie(isBoundary, alpha) = v->getBoundary();

    if ( isBoundary ) {

      // The stiffness matrix is a N by N matrix here
      // first the right_f_phi, move from left to right
      // then assign 0 to horizontal line and vertical line
      // except for (index, index)
      for ( size_t i = 0; i < index; ++ i ) {

        // right_f_phi, move K from left to right with alpha
        double K = _stiffness_matrix(index, i);
        _right_f_phi(i) -= K * alpha;

        // assign 0
        _stiffness_matrix (index, i) = 0;
        _stiffness_matrix (i, index) = 0;
      }

      for ( size_t i = index + 1; i < _stiffness_matrix.size1(); ++ i ) {
        double K = _stiffness_matrix(index, i);
        _right_f_phi(i) -= K * alpha;

        _stiffness_matrix (index, i) = 0;
        _stiffness_matrix (i, index) = 0;
      }

      // assign 1 to itself
      _stiffness_matrix (index, index) = 1;

      // assign alpha to itself
      _right_f_phi(index) = alpha;

      // Done
    }
  }

  FEM_triangular::FEM_triangular(std::pair<double, double> v1,
                                 std::pair<double, double> v2,
                                 size_t slice_x_num,
                                 size_t slice_y_num,
                                 utils::function p,
                                 utils::function q,
                                 utils::function f,
                                 size_t degree) :
    FEM_triangular(v1, v2, slice_x_num, slice_y_num, p, q, f,
                   degree, 2 /* this parameter is not important any more*/)
  { }

  std::pair<matrix<utils::function>, matrix<std::pair<size_t, size_t>>>
  FEM_triangular::K_element_matrix(const face &element_face) const {
    size_t size1, size2;
    size1 = size2 = _base_fun_vec.size();

    matrix<utils::function> N1(size1, size2);
    matrix<utils::function> N2(size1, size2);
    matrix<utils::function> N3(size1, size2);

    N1 = outer_prod(_base_diff_vec_kexi, _base_diff_vec_kexi);
    N2 = outer_prod(_base_diff_vec_kexi, _base_diff_vec_yita);
    N2 += outer_prod(_base_diff_vec_yita, _base_diff_vec_kexi);
    N3 = outer_prod(_base_diff_vec_yita, _base_diff_vec_yita);

    std::vector<size_t> vertices_index;
    // element x, y coordinates
    vector<double> x_coor, y_coor;

    std::tie(vertices_index, x_coor, y_coor) =
      __index_x_coor_y_coor_get(element_face);

    utils::function x_kexi, x_yita, y_kexi, y_yita;
    std::tie(x_kexi, x_yita, y_kexi, y_yita) =
      __local_derivative_function_get(x_coor, y_coor);

    utils::function detJe([=] (double x, double y) -> double {
                            auto result = fabs(x_kexi(x, y) * y_yita(x, y) -
                                               x_yita(x, y) * y_kexi(x, y));

                            if ( utils::almost_equal(result, (double)0) )
                              std::cout << "Got a zero here!\n";

                            return result;
                            // fabs(x_kexi(x, y) * y_yita(x, y) -
                            //      x_yita(x, y) * y_kexi(x, y));
                          });

    utils::function KY_to_x, KY_to_y;
    std::tie(KY_to_x, KY_to_y) = __linear_co_transform(x_coor, y_coor);

    utils::function g1e, g2e, g3e;

    g1e = [=] (double kexi, double yita) -> double {
            return
              _p(KY_to_x(kexi, yita),
                 KY_to_y(kexi, yita)) *
              (y_yita(kexi, yita) * y_yita(kexi, yita) +
               (x_yita(kexi, yita) * x_yita(kexi, yita)));
          };

    g2e = [=] (double kexi, double yita) -> double {
            return _p(KY_to_x(kexi, yita),
                      KY_to_y(kexi, yita)) *
              - (y_yita * y_kexi + x_yita * x_kexi)(kexi, yita);
          };

    g3e = [=] (double kexi, double yita) -> double {
            return _p(KY_to_x(kexi, yita),
                      KY_to_y(kexi, yita)) *
              (y_kexi * y_kexi + x_kexi * x_kexi)(kexi, yita);
          };

    auto function_matrix = (g1e * N1 + g2e * N2 + g3e * N3) / detJe;

    matrix<std::pair<size_t, size_t>> index_matrix(size1, size2);

    // assignment, index for assembly of global stiffness matrix
    for ( auto i = 0; i < size1; ++ i ) {
      for ( auto j = 0; j < size2; ++ j ) {
        // function
        //res.first(i, j) = tmp(i, j);
        // index
        std::pair<size_t, size_t> *tmp = &(index_matrix(i, j));
        tmp->first = vertices_index[i];
        tmp->second = vertices_index[j];
      }
    }

    return std::make_pair(std::move(function_matrix), std::move(index_matrix));

  }

  std::pair<matrix<utils::function>, matrix<std::pair<size_t, size_t>>>
  FEM_triangular::M_element_matrix(const face &element_face) const {
  size_t size1, size2;
    size1 = size2 = _base_fun_vec.size();

    auto res = std::make_pair(matrix<utils::function>(size1, size2),
                              matrix<std::pair<size_t, size_t>>(size1, size2));


    // element function matrix
    matrix<utils::function> NN(size1, size2);

    NN = outer_prod(_base_fun_vec, _base_fun_vec);

    std::vector<size_t> vertices_index;
    // element x, y coordinates
    vector<double> x_coor, y_coor;

    // get coordinate
    std::tie(vertices_index, x_coor, y_coor) =
      __index_x_coor_y_coor_get(element_face);

    // Get linear transformation function
    // Get differentiation functions
    utils::function KY_to_x, KY_to_y;
    utils::function x_kexi, x_yita, y_kexi, y_yita;
    std::tie(KY_to_x, KY_to_y) = __linear_co_transform(x_coor, y_coor);
    std::tie(x_kexi, x_yita, y_kexi, y_yita) =
      __local_derivative_function_get(x_coor, y_coor);

    utils::function q([=] (double kexi, double yita) -> double {
                        return _q(KY_to_x(kexi, yita),
                                  KY_to_y(kexi, yita));
                      });


    utils::function detJe([=] (double x, double y) -> double {
                            return
                              fabs(x_kexi(x, y) * y_yita(x, y) -
                                   x_yita(x, y) * y_kexi(x, y));
                          });

    auto tmp = q * NN * detJe;
    res.first = std::move(tmp);

    for ( size_t i = 0; i < size1; ++ i ) {
      for ( size_t j = 0; j < size2; ++ j ) {
        //res.first(i, j) = tmp(i, j);

        // index
        auto tmp = &(res.second(i, j));
        tmp->first = vertices_index[i];
        tmp->second = vertices_index[i];
      }
    }

    return res;
  }

  std::pair<vector<utils::function>, vector<size_t>>
  FEM_triangular::f_element_vector(const face &element_face) const {
    size_t size1 = _base_fun_vec.size();
    auto res = std::make_pair(vector<utils::function>(size1),
                              vector<size_t>(size1));

    auto N = _base_fun_vec;

    std::vector<size_t> vertices_index;
    vector<double> x_coor, y_coor;
    std::tie(vertices_index, x_coor, y_coor) =
      __index_x_coor_y_coor_get(element_face);

    utils::function KY_to_x, KY_to_y;
    std::tie(KY_to_x, KY_to_y) = __linear_co_transform(x_coor, y_coor);

    utils::function x_kexi, x_yita, y_kexi, y_yita;
    std::tie(x_kexi, x_yita, y_kexi, y_yita) =
      __local_derivative_function_get(x_coor, y_coor);

    utils::function f([=] (double kexi, double yita) -> double {
                        return _f(KY_to_x(kexi, yita),
                                  KY_to_y(kexi, yita));
                      });

    utils::function detJe([=] (double x, double y)-> double {
                            // std::cout << fabs(x_kexi(x, y) * y_yita(x, y) -
                            //                   x_yita(x, y) * y_kexi(x, y)) << '\n';
                            return fabs(x_kexi(x, y) * y_yita(x, y) -
                                        x_yita(x, y) * y_kexi(x, y));
                          });

    auto tmp = f * N * detJe;
    res.first = std::move(tmp);

    for ( size_t i = 0; i < size1; ++ i ) {
      //res.first(i) = tmp(i);

      auto tmp = &(res.second(i));
      *tmp = vertices_index[i];
    }

    return res;
  }

};
