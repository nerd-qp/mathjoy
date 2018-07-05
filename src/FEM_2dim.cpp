
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

  using boost::numeric::ublas::matrix;
  using boost::numeric::ublas::vector;

  // FEM_square solver's method
  std::pair<utils::function, utils::function>
  FEM_square::__linear_co_transform(const vector<double> &x_coor,
                                    const vector<double> &y_coor) const {
    utils::function KY_to_x, KY_to_y;

    // the base case
    // KY_to_x = KY_to_y = ([=] (double x, double y) -> double {
    //                        return 0;
    //                      });

    //std::cout << (_base_fun_vec[0] * x_coor[0]) (0, 1) << '\n';
    KY_to_x = utils::inner_prod(_base_fun_vec, x_coor);
    KY_to_y = utils::inner_prod(_base_fun_vec, y_coor);

    // for ( size_t i = 0; i < _base_fun_vec.size(); ++ i ) {
    //   KY_to_x += _base_fun_vec[i] * x_coor[i];
    //   KY_to_y += _base_fun_vec[i] * x_coor[i];
    // }

    //std::cout << (KY_to_y) (0, 1) << "here\n";
    return std::make_pair(KY_to_x, KY_to_y);
  }

  inline std::tuple<utils::function,
                    utils::function,
                    utils::function,
                    utils::function>
  FEM_square::__local_derivative_function_get(const vector<double> &x_coor,
                                              const vector<double> &y_coor) const {
    const auto &N_kexi = _base_diff_vec_kexi;
    const auto &N_yita = _base_diff_vec_yita;

    // Get the points out, and make derivative functions
    utils::function x_kexi, x_yita, y_kexi, y_yita;

    // for ( auto i : x_coor )
    //   std::cout << i << '\n';

    x_kexi = utils::inner_prod(N_kexi, x_coor);
    x_yita = utils::inner_prod(N_yita, x_coor);
    y_kexi = utils::inner_prod(N_kexi, y_coor);
    y_yita = utils::inner_prod(N_yita, y_coor);

    // std::cout << x_yita(-1, 0) << "x_kexi\n";

    return std::make_tuple(x_kexi, x_yita,
                           y_kexi, y_yita);
  }


  std::tuple<std::vector<size_t>, vector<double>, vector<double>>
  FEM_square::__index_x_coor_y_coor_get(const face &element_face) const {
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

  void
  FEM_square::__assemble_local_stiffness_matrix(const face &element_face) {
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

  void
  FEM_square::__assemble_local_right_f_phi(const face &element_face) {
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

  // Delegating constructors
  FEM_square::FEM_square (std::pair<double, double> v1,
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
    _base_diff_vec_yita(freedom)
  { }

  void
  FEM_square::assemble_stiffness_matrix() {
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

    //std::cout << '\n';

  }

  void
  FEM_square::assemble_right_f_phi() {
    auto all_faces = _mesh->get_faces();

    boost::progress_display show_progress(all_faces.size());
    for ( size_t i = 0; i < all_faces.size(); ++ i ) {
      __assemble_local_right_f_phi(all_faces[i]);

      // have to use flush
      // if ( (i * 100 / all_faces.size()) % 10 == 0 )
      //   std::cout << i << '.' << std::flush;

      ++show_progress;
    }

    // for ( const auto &face : all_faces ) {

    //   std::cout << '.';
    // }

    // std::cout << '\n';

  }

  void
  FEM_square::handle_boundary_point(size_t index) {
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

  FEM_square::FEM_square (std::pair<double, double> v1,
                          std::pair<double, double> v2,
                          size_t slice_x_num,
                          size_t slice_y_num,
                          utils::function p,
                          utils::function q,
                          utils::function f,
                          size_t degree) :
    FEM_square(v1, v2, slice_x_num, slice_y_num, p, q, f,
               degree, (pow(degree + 1, 2)))
  {
    if ( degree > 1 )
      throw "only support degree of 1 now!\n";
    //_mesh = std::unique_ptr<simple_mesh>(new simple_square_mesh((v1, v2, slice_x_num, slice_y_num)));
    _mesh = utils::make_unique<simple_square_mesh>(v1, v2,
                                                   slice_x_num,
                                                   slice_y_num);
    _quadrature_table =
      utils::make_unique<two_dimension::canonical::square>(5);

    //_base_fun_vec = Lagrange_gen(degree);

    /*
     * NOTE evil code block here...
     */
    double range_a = -1, range_b = 1;
    // generate distance vector
    std::vector<double> distance(degree,
                                 (range_b - range_a) / (double)degree);

    // generate interpolation nodes, would throw exception
    auto interpolation_nodes = utils::interpolation_gen(range_a,
                                                        range_b,
                                                        distance);

    // std::vector<double> interpolation_nodes_x, interpolation_nodes_y;
    // interpolation_nodes_x.reserve(distance.size() * distance.size());
    // interpolation_nodes_y.reserve(distance.size() * distance.size());

    // utils::fully_connect(tmp, tmp, [&](double x, double y) {
    //                           interpolation_nodes_x.emplace_back(x);
    //                           interpolation_nodes_y.emplace_back(y);
    //                         });

    // utils::zip(interpolation_nodes_x.begin(), interpolation_nodes_x.end(),
    //     interpolation_nodes_y.begin(), interpolation_nodes_y.end(),
    //     [] (double x, double y) {
    //       std::cout << x << ", " << y << '\n';
    //     });

    // std::cout << '\n';
    // generate 1 dim LagrangePolynomial and its corresponding differentiation
    // auto tmp1 = utils::LagrangePolynomial(interpolation_nodes);
    // auto tmp2 = utils::LagrangePolynomial_diff(interpolation_nodes);

    // vector<utils::function_1dim> onedim_functions(tmp1.size());
    // vector<utils::function_1dim> onedim_diff(tmp2.size());

    // for ( size_t i = 0; i < onedim_functions.size(); ++ i ) {
    //   onedim_functions(i) = tmp1[i];
    //   onedim_diff(i) = tmp2[i];
    // }

    // _base_diff_vec_kexi = utils::tensor_product_fun(onedim_diff, onedim_functions);
    // _base_diff_vec_yita = utils::tensor_product_fun(onedim_functions, onedim_diff);

    // _base_fun_vec =
    //   utils::tensor_product_fun(onedim_functions, onedim_functions);

    /*
     * NOTE Evil code block, should be a function here
     */

    std::tie(_base_fun_vec, _base_diff_vec_kexi, _base_diff_vec_yita) =
      utils::Lagrange_with_diff_2dim(interpolation_nodes, interpolation_nodes);
    // auto res = utils::Lagrange_with_diff_2dim(interpolation_nodes, interpolation_nodes);

    // std::cout << "old : " << _base_fun_vec[0](10, 10) << std::endl;
    // _base_fun_vec = std::get<0>(res);

    assemble_stiffness_matrix();
    assemble_right_f_phi();

    for ( size_t i = 0; i < _mesh->get_vertices().size(); ++ i )
      handle_boundary_point(i);

    // auto tmp1 = K_element_matrix(_mesh->get_faces()[0]);
    // auto tmp2 = M_element_matrix(_mesh->get_faces()[0]);
    // auto tmp3 = f_element_vector(_mesh->get_faces()[0]);

    // std::cout << (_quadrature_table->quadrature_compute(tmp1.first)) << '\n'
    //           << (_quadrature_table->quadrature_compute(tmp2.first)) << '\n'
    //           << (_quadrature_table->quadrature_compute(tmp3.first)) << '\n'
    //           << "\nVery exciting isn't it??\n";
  }

  std::pair<matrix<utils::function>, matrix<std::pair<size_t, size_t>>>
  FEM_square::K_element_matrix(const face &element_face) const {
    // Initilize with sizes
    size_t size1, size2;
    size1 = size2 = _base_fun_vec.size();

    // element function matrix
    matrix<utils::function> N1(size1, size2);
    matrix<utils::function> N2(size1, size2);
    matrix<utils::function> N3(size1, size2);

    // std::cout << _base_diff_vec_kexi(0)(1, 1) << '\n'
    //           << _base_diff_vec_kexi(1)(1, 1) << '\n'
    //           << _base_diff_vec_kexi(2)(1, 1) << '\n'
    //           << _base_diff_vec_kexi(3)(1, 1) << '\n';

    N1 = outer_prod(_base_diff_vec_kexi, _base_diff_vec_kexi);
    N2 = outer_prod(_base_diff_vec_kexi, _base_diff_vec_yita);
    N2 += outer_prod(_base_diff_vec_yita, _base_diff_vec_kexi);
    N3 = outer_prod(_base_diff_vec_yita, _base_diff_vec_yita);

    std::vector<size_t> vertices_index;
    // element x, y coordinates
    vector<double> x_coor, y_coor;

    std::tie(vertices_index, x_coor, y_coor) =
      __index_x_coor_y_coor_get(element_face);

    // for ( size_t i = 0; i < x_coor.size(); ++ i ) {
    //   std::cout << x_coor[i] << ", " << y_coor[i] << '\n';
    // }
    // std::cout << "========================================\n";

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

    //std::cout << detJe(-1, -1) << "detJe\n";

    // Kexi Yita to X, x(kexi, yita)
    // Kexi Yita to Y, y(kexi, yita)
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

    // Some optimization to lower operations (decrypted), useless
    // auto g1 = g1e / detJe;
    // auto g2 = g2e / detJe;
    // auto g3 = g3e / detJe;

    // auto tmp = (g1 * N1 + g2 * N2 + g3 * N3);
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

    return std::make_pair(function_matrix, index_matrix);
  }

  std::pair<matrix<utils::function>, matrix<std::pair<size_t, size_t>>>
  FEM_square::M_element_matrix(const face &element_face) const {
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
  FEM_square::f_element_vector(const face &element_face) const {
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

}
