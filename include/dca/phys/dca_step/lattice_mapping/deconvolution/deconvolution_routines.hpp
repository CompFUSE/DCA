// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides routines for the deconvolution step.

#ifndef DCA_PHYS_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_DECONVOLUTION_ROUTINES_HPP
#define DCA_PHYS_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_DECONVOLUTION_ROUTINES_HPP

#include <cmath>
#include <complex>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/math/function_transform/basis_transform/basis_transform.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_sp.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"

namespace dca {
namespace phys {
namespace latticemapping {
// dca::phys::latticemapping::

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
class deconvolution_routines {
public:
  using target_k_cluster_type = typename target_k_dmn_t::parameter_type;
  using target_r_cluster_type = typename target_k_cluster_type::dual_type;
  using target_r_dmn_t = func::dmn_0<target_r_cluster_type>;
  using trafo_k_to_r_type = math::transform::basis_transform<target_k_dmn_t, target_r_dmn_t>;
  using trafo_r_to_k_type = math::transform::basis_transform<target_r_dmn_t, target_k_dmn_t>;

public:
  deconvolution_routines(parameters_type& parameters_ref);

  void compute_T_inv_matrix(double epsilon, dca::linalg::Matrix<double, dca::linalg::CPU>& T_eps);
  void compute_T_inv_matrix(double epsilon,
                            dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU>& T_eps);

  const auto& get_T() const {
    return T_;
  }
  const auto& get_T_symmetrized() const {
    return T_symmetrized_;
  }

  const auto& get_T_source() const {
    return T_source_;
  }
  const auto& get_T_source_symmetrized() const {
    return T_source_symmetrized_;
  }

private:
  void compute_phi_inv(double epsilon);

  template <typename LhsKDomain>
  void initializeProjectionOperator(const func::function<double, target_r_dmn_t>& phi_r,
                                    linalg::Matrix<double, dca::linalg::CPU>& projection_op);

private:
  parameters_type& parameters;

protected:
  func::function<double, target_r_dmn_t> phi_r_inv;

private:
  func::function<double, target_r_dmn_t> phi_r_;

  dca::linalg::Matrix<double, dca::linalg::CPU> T_;
  dca::linalg::Matrix<double, dca::linalg::CPU> T_symmetrized_;

  dca::linalg::Matrix<double, dca::linalg::CPU> T_source_;
  dca::linalg::Matrix<double, dca::linalg::CPU> T_source_symmetrized_;
};

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::deconvolution_routines(
    parameters_type& parameters_ref)
    : parameters(parameters_ref),

      phi_r_inv("phi_r_inv"),
      phi_r_("phi(r)"),

      T_("T", std::make_pair(target_k_dmn_t::dmn_size(), target_k_dmn_t::dmn_size())),
      T_symmetrized_("T-symmetrized",
                     std::make_pair(target_k_dmn_t::dmn_size(), target_k_dmn_t::dmn_size())),

      T_source_("T-source", std::make_pair(source_k_dmn_t::dmn_size(), target_k_dmn_t::dmn_size())),
      T_source_symmetrized_("T-source-symmetrized",
                            std::make_pair(source_k_dmn_t::dmn_size(), target_k_dmn_t::dmn_size())) {
  clustermapping::CoarsegrainingSp<parameters_type> coarsegrain_obj(parameters);

  coarsegrain_obj.compute_phi_r(phi_r_);

  func::function<double, target_r_dmn_t> phi_r_symmetrized(phi_r_, "phi_symmetrized(r)");
  symmetrize::execute(phi_r_symmetrized);

  // Compute target (lattice) k-domain to target k-domain projection operators.
  initializeProjectionOperator<target_k_dmn_t>(phi_r_, T_);
  initializeProjectionOperator<target_k_dmn_t>(phi_r_symmetrized, T_symmetrized_);

  // Compute target (lattice) k-domain to source (cluster) k-domain projection operators.
  initializeProjectionOperator<source_k_dmn_t>(phi_r_, T_source_);
  initializeProjectionOperator<source_k_dmn_t>(phi_r_symmetrized, T_source_symmetrized_);
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
template <typename LhsKDomain>
void deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::initializeProjectionOperator(
    const func::function<double, target_r_dmn_t>& phi_r, linalg::Matrix<double, dca::linalg::CPU>& T) {
  using trafo_r_to_lhs_k_type = math::transform::basis_transform<target_r_dmn_t, LhsKDomain>;

  const linalg::Matrix<std::complex<double>, dca::linalg::CPU>& T_k_to_r =
      trafo_k_to_r_type::get_transformation_matrix();
  const linalg::Matrix<std::complex<double>, dca::linalg::CPU>& T_r_to_k =
      trafo_r_to_lhs_k_type::get_transformation_matrix();

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> T_k_to_r_scaled(T_k_to_r,
                                                                              "T_k_to_r_scaled");
  for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
    for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
      T_k_to_r_scaled(i, j) *= phi_r(i);

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> T_k_to_k(
      "T_k_to_k", std::make_pair(LhsKDomain::dmn_size(), target_k_dmn_t::dmn_size()));

  linalg::matrixop::gemm(T_r_to_k, T_k_to_r_scaled, T_k_to_k);

  for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
    for (int i = 0; i < LhsKDomain::dmn_size(); i++)
      T(i, j) = std::real(T_k_to_k(i, j));

  // Normalize the rows.
  for (int i = 0; i < LhsKDomain::dmn_size(); i++) {
    double result = 0;

    for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
      result += T(i, j);

    for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
      T(i, j) /= result;
  }
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
void deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::compute_phi_inv(
    double epsilon) {
  for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
    phi_r_inv(i) = std::abs(phi_r_(i)) > epsilon ? 1. / phi_r_(i) : 0.;  // 1./epsilon;
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
void deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::compute_T_inv_matrix(
    double epsilon, dca::linalg::Matrix<double, dca::linalg::CPU>& T_eps) {
  compute_phi_inv(epsilon);

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU>& T_k_to_r =
      trafo_k_to_r_type::get_transformation_matrix();
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU>& T_r_to_k =
      trafo_r_to_k_type::get_transformation_matrix();

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> T_k_to_r_scaled("T_k_to_r_scaled");
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> T_k_to_k("T_k_to_r");

  T_k_to_k.resize(T_k_to_r.size());
  T_k_to_r_scaled = T_k_to_r;

  for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
    for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
      T_k_to_r_scaled(i, j) *= phi_r_inv(i);

  dca::linalg::matrixop::gemm(T_r_to_k, T_k_to_r_scaled, T_k_to_k);

  for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
    for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
      T_eps(i, j) = real(T_k_to_k(i, j));
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
void deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::compute_T_inv_matrix(
    double epsilon, dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU>& T_eps) {
  compute_phi_inv(epsilon);

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU>& T_k_to_r =
      trafo_k_to_r_type::get_transformation_matrix();
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU>& T_r_to_k =
      trafo_r_to_k_type::get_transformation_matrix();

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> T_k_to_r_scaled("T_k_to_r_scaled");
  T_k_to_r_scaled = T_k_to_r;

  for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
    for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
      T_k_to_r_scaled(i, j) *= phi_r_inv(i);

  dca::linalg::matrixop::gemm(T_r_to_k, T_k_to_r_scaled, T_eps);
}

}  // latticemapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_DECONVOLUTION_ROUTINES_HPP
