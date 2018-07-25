// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the interpolation matrices for the coarsegraining.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_INTERPOLATION_MATRICES_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_INTERPOLATION_MATRICES_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <mutex>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/math/function_transform/basis_transform/basis_transform.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegrain_domain_names.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_domain.hpp"
#include "dca/phys/domains/cluster/centered_cluster_domain.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename scalar_type, typename k_dmn, typename K_dmn>
class interpolation_matrices {};

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
class interpolation_matrices<scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>> {
public:
  using r_dmn = typename k_dmn::parameter_type::dual_type;

  using q_dmn = func::dmn_0<coarsegraining_domain<K_dmn, NAME>>;
  using r_centered_dmn = func::dmn_0<domains::centered_cluster_domain<r_dmn>>;

  using trafo_k_to_r_type = math::transform::basis_transform<k_dmn, r_centered_dmn>;
  using trafo_r_to_q_type = math::transform::basis_transform<r_centered_dmn, q_dmn>;
  using trafo_matrix_type = typename trafo_k_to_r_type::matrix_type;

  using matrix_type = dca::linalg::Matrix<scalar_type, dca::linalg::CPU>;

public:
  static func::function<matrix_type, K_dmn>& get() {
    static func::function<matrix_type, K_dmn> k_to_q("k_to_q (" +
                                                     q_dmn::parameter_type::get_name() + ")");
    return k_to_q;
  }

  static matrix_type& get(int k_ind) {
    static func::function<matrix_type, K_dmn>& k_to_q = get();
    return k_to_q(k_ind);
  }

  static bool is_initialized() {
    return initialized_;
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

  static bool initialized_;
  static std::mutex initialization_mutex_;
};
template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
bool interpolation_matrices<scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>>::initialized_ =
    false;
template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
std::mutex interpolation_matrices<
    scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>>::initialization_mutex_;

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename concurrency_type>
void interpolation_matrices<scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>>::resize_matrices(
    concurrency_type& concurrency) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t interpolation-matrices " << to_str(NAME) << " initialization started ... ";

  for (int K_ind = 0; K_ind < K_dmn::dmn_size(); K_ind++) {
    matrix_type& T_k_to_q = get(K_ind);

    T_k_to_q.resizeNoCopy(std::pair<int, int>(q_dmn::dmn_size(), k_dmn::dmn_size()));

    for (int j = 0; j < k_dmn::dmn_size(); j++)
      for (int i = 0; i < q_dmn::dmn_size(); i++)
        T_k_to_q(i, j) = 0;
  }

  r_centered_dmn::parameter_type::initialize();
}

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename concurrency_type>
void interpolation_matrices<scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>>::print_memory_used(
    concurrency_type& concurrency) {
  if (concurrency.id() == concurrency.first()) {
    std::pair<int, int> capacity = get(0).capacity();
    std::cout << " stopped ( "
              << sizeof(scalar_type) * capacity.first * capacity.second * 1.e-6 * K_dmn::dmn_size()
              << " Mbytes) \n\n";
  }
}

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename concurrency_type>
void interpolation_matrices<scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>>::initialize(
    concurrency_type& concurrency) {
  assert(NAME == K or NAME == TETRAHEDRON_K);

  std::unique_lock<std::mutex> lock(initialization_mutex_);

  if (is_initialized())
    return;

  resize_matrices(concurrency);

  K_dmn K_dmn_obj;
  std::pair<int, int> bounds = concurrency.get_bounds(K_dmn_obj);

  trafo_matrix_type trafo_k_to_q;
  trafo_k_to_q.resizeNoCopy(std::pair<int, int>(q_dmn::dmn_size(), k_dmn::dmn_size()));

  for (int K_ind = bounds.first; K_ind < bounds.second; K_ind++) {
    {
      q_dmn::parameter_type::set_elements(K_ind);

      // trafo_k_to_r_type::is_initialized() = false;
      trafo_r_to_q_type::is_initialized() = false;

      trafo_matrix_type& trafo_r_to_q = trafo_r_to_q_type::get_transformation_matrix();
      trafo_matrix_type& trafo_k_to_r = trafo_k_to_r_type::get_transformation_matrix();

      dca::linalg::matrixop::gemm(trafo_r_to_q, trafo_k_to_r, trafo_k_to_q);
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

  initialized_ = true;
  print_memory_used(concurrency);
}

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename concurrency_type>
void interpolation_matrices<scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>>::initialize(
    concurrency_type& concurrency, int Q_ind) {
  assert(NAME == K_PLUS_Q or NAME == Q_MINUS_K);

  std::unique_lock<std::mutex> lock(initialization_mutex_);

  if (is_initialized())
    return;

  resize_matrices(concurrency);

  K_dmn K_dmn_obj;
  std::pair<int, int> bounds = concurrency.get_bounds(K_dmn_obj);

  trafo_matrix_type trafo_k_to_q;
  trafo_k_to_q.resizeNoCopy(std::pair<int, int>(q_dmn::dmn_size(), k_dmn::dmn_size()));

  for (int K_ind = bounds.first; K_ind < bounds.second; K_ind++) {
    {
      q_dmn::parameter_type::set_elements(K_ind, Q_ind);

      // trafo_k_to_r_type::is_initialized() = false;
      trafo_r_to_q_type::is_initialized() = false;

      trafo_matrix_type& trafo_r_to_q = trafo_r_to_q_type::get_transformation_matrix();
      trafo_matrix_type& trafo_k_to_r = trafo_k_to_r_type::get_transformation_matrix();

      dca::linalg::matrixop::gemm(trafo_r_to_q, trafo_k_to_r, trafo_k_to_q);
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

  initialized_ = true;
  print_memory_used(concurrency);
}

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename scalar_type_1, typename scalar_type_2>
void interpolation_matrices<scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>>::cast(
    scalar_type_1& x, scalar_type_2& y) {
  x = y;
}

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
template <typename scalar_type_1, typename scalar_type_2>
void interpolation_matrices<scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>>::cast(
    scalar_type_1& x, std::complex<scalar_type_2>& y) {
  assert(std::abs(std::imag(y)) < 1.e-6);
  x = std::real(y);
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_INTERPOLATION_MATRICES_HPP
