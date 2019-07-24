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

#include "dca/config/profiler.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrix_view.hpp"
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

  using Matrix = dca::linalg::Matrix<scalar_type, dca::linalg::CPU>;
  using MatrixView = dca::linalg::MatrixView<scalar_type, dca::linalg::CPU>;

  static MatrixView get(int k_ind) {
    auto& k_to_q = get();
    return MatrixView(&k_to_q(0, 0, k_ind), std::make_pair(q_dmn::dmn_size(), k_dmn::dmn_size()));
  }

  static void set_q_idx(int q_idx) {
    assert(q_idx >= 0 && q_idx < q_dmn::dmn_size());
    q_idx_ = q_idx;
  }

private:
  // Interpolation matrices concatenated column-wise.
  using Matrices = func::function<scalar_type, func::dmn_variadic<q_dmn, k_dmn, K_dmn>>;
  static auto get() -> Matrices&;

  static auto initialize() -> Matrices;

  template <typename scalar_type_1, typename scalar_type_2>
  inline static void cast(scalar_type_1& x, scalar_type_2& y);

  template <typename scalar_type_1, typename scalar_type_2>
  inline static void cast(scalar_type_1& x, std::complex<scalar_type_2>& y);

  static int q_idx_;
};
template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
int interpolation_matrices<scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>>::q_idx_ =
    -1;

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
auto interpolation_matrices<scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>>::get()
    -> Matrices& {
  static Matrices k_to_q(initialize(), "k_to_q (" + q_dmn::parameter_type::get_name() + ")");
  return k_to_q;
}

template <typename scalar_type, typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
auto interpolation_matrices<scalar_type, k_dmn, func::dmn_0<coarsegraining_domain<K_dmn, NAME>>>::initialize()
    -> Matrices {
  Profiler profiler(__FUNCTION__, "Interpolation matrices", __LINE__);

  Matrices matrices;
  r_centered_dmn::parameter_type::initialize();

  trafo_matrix_type trafo_k_to_q;
  trafo_k_to_q.resizeNoCopy(std::pair<int, int>(q_dmn::dmn_size(), k_dmn::dmn_size()));

  for (int K_ind = 0; K_ind < K_dmn::dmn_size(); ++K_ind) {
    if (NAME == K || NAME == TETRAHEDRON_K) {
      q_dmn::parameter_type::set_elements(K_ind);
    }
    else if (NAME == K_PLUS_Q || NAME == Q_MINUS_K) {
      if (q_idx_ == -1)
        throw(std::logic_error("q index not set."));
      q_dmn::parameter_type::set_elements(K_ind, q_idx_);
    }
    else {
      throw(std::logic_error("Non matching NAME."));
    }

    trafo_r_to_q_type::is_initialized() = false;

    trafo_matrix_type& trafo_r_to_q = trafo_r_to_q_type::get_transformation_matrix();
    trafo_matrix_type& trafo_k_to_r = trafo_k_to_r_type::get_transformation_matrix();

    dca::linalg::matrixop::gemm(trafo_r_to_q, trafo_k_to_r, trafo_k_to_q);

    MatrixView T_k_to_q =
        MatrixView(&matrices(0, 0, K_ind), std::make_pair(q_dmn::dmn_size(), k_dmn::dmn_size()));

    for (int j = 0; j < k_dmn::dmn_size(); j++)
      for (int i = 0; i < q_dmn::dmn_size(); i++)
        cast(T_k_to_q(i, j), trafo_k_to_q(i, j));
  }

  return matrices;
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

}  // namespace clustermapping
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_INTERPOLATION_MATRICES_HPP
