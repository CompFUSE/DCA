//-*-C++-*-

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_INTERPOLATION_MATRICES_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_INTERPOLATION_MATRICES_H

namespace DCA {
template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
class interpolation_matrices<scalar_type, k_dmn, dmn_0<coarsegraining_domain<K_dmn, NAME>>> {
  typedef typename k_dmn::parameter_type::dual_type r_dmn;

  typedef dmn_0<coarsegraining_domain<K_dmn, NAME>> q_dmn;
  typedef dmn_0<centered_cluster_domain<r_dmn>> r_centered_dmn;

  typedef math_algorithms::functional_transforms::basis_transform<k_dmn, r_centered_dmn> trafo_k_to_r_type;
  typedef math_algorithms::functional_transforms::basis_transform<r_centered_dmn, q_dmn> trafo_r_to_q_type;

  typedef typename trafo_k_to_r_type::matrix_type trafo_matrix_type;

public:
  typedef LIN_ALG::matrix<scalar_type, LIN_ALG::CPU> matrix_type;

public:
  static FUNC_LIB::function<matrix_type, K_dmn>& get() {
    assert(is_initialized() == true);

    static FUNC_LIB::function<matrix_type, K_dmn> k_to_q("k_to_q (" +
                                                         q_dmn::parameter_type::get_name() + ")");
    return k_to_q;
  }

  static matrix_type& get(int k_ind) {
    assert(is_initialized() == true);

    static FUNC_LIB::function<matrix_type, K_dmn>& k_to_q = get();
    return k_to_q(k_ind);
  }

  static bool& is_initialized() {
    static bool initialized = false;
    return initialized;
  }

  template <typename concurrency_type>
  static void initialize(concurrency_type& concurrency);

  template <typename concurrency_type>
  static void initialize(concurrency_type& concurrency, int Q_ind);

private:
  template <typename concurrency_type>
  static void resize_matrices(concurrency_type& concurrency);

  template <typename concurrency_type>
  static void print_memory_used(concurrency_type& concurrency);

  template <typename scalar_type_1, typename scalar_type_2>
  inline static void cast(scalar_type_1& x, scalar_type_2& y);

  template <typename scalar_type_1, typename scalar_type_2>
  inline static void cast(scalar_type_1& x, std::complex<scalar_type_2>& y);
};

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename concurrency_type>
void interpolation_matrices<scalar_type, k_dmn, dmn_0<coarsegraining_domain<K_dmn, NAME>>>::resize_matrices(
    concurrency_type& concurrency) {
  if (concurrency.id() == 0)
    std::cout << "\n\n\t interpolation-matrices " << to_str(NAME) << " initialization started ... ";

  is_initialized() = true;

  for (int K_ind = 0; K_ind < K_dmn::dmn_size(); K_ind++) {
    matrix_type& T_k_to_q = get(K_ind);

    T_k_to_q.resize_no_copy(std::pair<int, int>(q_dmn::dmn_size(), k_dmn::dmn_size()));

    for (int j = 0; j < k_dmn::dmn_size(); j++)
      for (int i = 0; i < q_dmn::dmn_size(); i++)
        T_k_to_q(i, j) = 0;
  }

  r_centered_dmn::parameter_type::initialize();
}

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename concurrency_type>
void interpolation_matrices<scalar_type, k_dmn, dmn_0<coarsegraining_domain<K_dmn, NAME>>>::print_memory_used(
    concurrency_type& concurrency) {
  if (concurrency.id() == 0) {
    std::pair<int, int> global_size = get(0).get_global_size();
    std::cout << " stopped ( "
              << sizeof(scalar_type) * global_size.first * global_size.second * 1.e-6 *
                     K_dmn::dmn_size()
              << " Mbytes) \n\n";
  }
}

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename concurrency_type>
void interpolation_matrices<scalar_type, k_dmn, dmn_0<coarsegraining_domain<K_dmn, NAME>>>::initialize(
    concurrency_type& concurrency) {
  assert(NAME == K or NAME == TETRAHEDRON_K);

  if (is_initialized())
    return;

  resize_matrices(concurrency);

  K_dmn K_dmn_obj;
  std::pair<int, int> bounds = concurrency.get_bounds(K_dmn_obj);

  trafo_matrix_type trafo_k_to_q;
  trafo_k_to_q.resize_no_copy(std::pair<int, int>(q_dmn::dmn_size(), k_dmn::dmn_size()));

  for (int K_ind = bounds.first; K_ind < bounds.second; K_ind++) {
    {
      q_dmn::parameter_type::set_elements(K_ind);

      // trafo_k_to_r_type::is_initialized() = false;
      trafo_r_to_q_type::is_initialized() = false;

      trafo_matrix_type& trafo_r_to_q = trafo_r_to_q_type::get_transformation_matrix();
      trafo_matrix_type& trafo_k_to_r = trafo_k_to_r_type::get_transformation_matrix();

      for (int j = 0; j < r_centered_dmn::dmn_size(); j++)
        for (int i = 0; i < q_dmn::dmn_size(); i++)
          trafo_r_to_q(i, j) *= r_centered_dmn::parameter_type::get_weights()[j];

      LIN_ALG::GEMM<LIN_ALG::CPU>::execute(trafo_r_to_q, trafo_k_to_r, trafo_k_to_q);
    }

    {
      matrix_type& T_k_to_q = get(K_ind);

      T_k_to_q.resize(std::pair<int, int>(q_dmn::dmn_size(), k_dmn::dmn_size()));

      for (int j = 0; j < k_dmn::dmn_size(); j++)
        for (int i = 0; i < q_dmn::dmn_size(); i++)
          cast(T_k_to_q(i, j), trafo_k_to_q(i, j));
    }
  }

  for (int K_ind = 0; K_ind < K_dmn::dmn_size(); K_ind++)
    concurrency.sum(get(K_ind));

  print_memory_used(concurrency);
}

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename concurrency_type>
void interpolation_matrices<scalar_type, k_dmn, dmn_0<coarsegraining_domain<K_dmn, NAME>>>::initialize(
    concurrency_type& concurrency, int Q_ind) {
  assert(NAME == K_PLUS_Q or NAME == Q_MINUS_K);

  if (is_initialized())
    return;

  resize_matrices(concurrency);

  K_dmn K_dmn_obj;
  std::pair<int, int> bounds = concurrency.get_bounds(K_dmn_obj);

  trafo_matrix_type trafo_k_to_q;
  trafo_k_to_q.resize_no_copy(std::pair<int, int>(q_dmn::dmn_size(), k_dmn::dmn_size()));

  for (int K_ind = bounds.first; K_ind < bounds.second; K_ind++) {
    {
      q_dmn::parameter_type::set_elements(K_ind, Q_ind);

      // trafo_k_to_r_type::is_initialized() = false;
      trafo_r_to_q_type::is_initialized() = false;

      trafo_matrix_type& trafo_r_to_q = trafo_r_to_q_type::get_transformation_matrix();
      trafo_matrix_type& trafo_k_to_r = trafo_k_to_r_type::get_transformation_matrix();

      for (int j = 0; j < r_centered_dmn::dmn_size(); j++)
        for (int i = 0; i < q_dmn::dmn_size(); i++)
          trafo_r_to_q(i, j) *= r_centered_dmn::parameter_type::get_weights()[j];

      LIN_ALG::GEMM<LIN_ALG::CPU>::execute(trafo_r_to_q, trafo_k_to_r, trafo_k_to_q);
    }

    {
      matrix_type& T_k_to_q = get(K_ind);

      T_k_to_q.resize(std::pair<int, int>(q_dmn::dmn_size(), k_dmn::dmn_size()));

      for (int j = 0; j < k_dmn::dmn_size(); j++)
        for (int i = 0; i < q_dmn::dmn_size(); i++)
          cast(T_k_to_q(i, j), trafo_k_to_q(i, j));
    }
  }

  for (int K_ind = 0; K_ind < K_dmn::dmn_size(); K_ind++)
    concurrency.sum(get(K_ind));

  print_memory_used(concurrency);
}

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename scalar_type_1, typename scalar_type_2>
void interpolation_matrices<scalar_type, k_dmn, dmn_0<coarsegraining_domain<K_dmn, NAME>>>::cast(
    scalar_type_1& x, scalar_type_2& y) {
  x = y;
}

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename scalar_type_1, typename scalar_type_2>
void interpolation_matrices<scalar_type, k_dmn, dmn_0<coarsegraining_domain<K_dmn, NAME>>>::cast(
    scalar_type_1& x, std::complex<scalar_type_2>& y) {
  assert(std::abs(imag(y)) < 1.e-6);
  x = real(y);
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_INTERPOLATION_MATRICES_H
