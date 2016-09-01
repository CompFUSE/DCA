// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class provides routines for the deconvolution step.

#ifndef PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_ROUTINES_H
#define PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_ROUTINES_H

#include <cmath>
#include <complex>
#include <utility>

#include "comp_library/function_library/include_function_library.h"
#include "comp_library/linalg/linalg.hpp"
#include "math_library/functional_transforms/basis_transforms/basis_transform.h"
#include "phys_library/DCA+_step/cluster_mapping/coarsegraining_step/coarsegraining_sp.h"
#include "phys_library/DCA+_step/symmetrization/symmetrize.h"

namespace DCA {

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
class deconvolution_routines {
public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using target_k_cluster_type = typename target_k_dmn_t::parameter_type;
  using target_r_cluster_type = typename target_k_cluster_type::dual_type;
  using target_r_dmn_t = dmn_0<target_r_cluster_type>;
  using trafo_k_to_r_type =
      math_algorithms::functional_transforms::basis_transform<target_k_dmn_t, target_r_dmn_t>;
  using trafo_r_to_k_type =
      math_algorithms::functional_transforms::basis_transform<target_r_dmn_t, target_k_dmn_t>;

public:
  deconvolution_routines(parameters_type& parameters_ref);

  void compute_T_inv_matrix(double epsilon, dca::linalg::Matrix<double, dca::linalg::CPU>& T_eps);
  void compute_T_inv_matrix(double epsilon,
                            dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU>& T_eps);

private:
  void initialize();

  void compute_phi_inv(double epsilon);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

protected:
  FUNC_LIB::function<double, target_r_dmn_t> phi_r;
  FUNC_LIB::function<double, target_r_dmn_t> phi_r_symmetrized;

  FUNC_LIB::function<double, target_r_dmn_t> phi_r_inv;

  dca::linalg::Matrix<double, dca::linalg::CPU> T;
  dca::linalg::Matrix<double, dca::linalg::CPU> T_symmetrized;
};

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::deconvolution_routines(
    parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      phi_r("phi(r)"),
      phi_r_symmetrized("phi_{sym}(r)"),

      phi_r_inv("phi_r_inv"),

      T("T            (deconvolution_routines)",
        std::pair<int, int>(target_k_dmn_t::dmn_size(), target_k_dmn_t::dmn_size())),
      T_symmetrized("T_symmetrize (deconvolution_routines)",
                    std::pair<int, int>(target_k_dmn_t::dmn_size(), target_k_dmn_t::dmn_size())) {
  initialize();
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
void deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::initialize() {
  DCA::coarsegraining_sp<parameters_type, source_k_dmn_t> coarsegrain_obj(parameters);

  coarsegrain_obj.compute_phi_r(phi_r);

  phi_r_symmetrized = phi_r;

  symmetrize::execute(phi_r_symmetrized);

  {
    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU>& T_k_to_r =
        trafo_k_to_r_type::get_transformation_matrix();
    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU>& T_r_to_k =
        trafo_r_to_k_type::get_transformation_matrix();

    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> T_k_to_r_scaled("T_k_to_r_scaled");
    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> T_k_to_k("T_k_to_r");

    T_k_to_k.copy_from(T_k_to_r);  // resize the matrix;

    {
      T_k_to_r_scaled.copy_from(T_k_to_r);

      for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
        for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
          T_k_to_r_scaled(i, j) *= phi_r(i);

      dca::linalg::matrixop::gemm(T_r_to_k, T_k_to_r_scaled, T_k_to_k);

      for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
        for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
          T(i, j) = real(T_k_to_k(i, j));

      for (int i = 0; i < target_k_dmn_t::dmn_size(); i++) {
        double result = 0;

        for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
          result += T(i, j);

        for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
          T(i, j) /= result;
      }
    }

    {
      T_k_to_r_scaled.copy_from(T_k_to_r);

      for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
        for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
          T_k_to_r_scaled(i, j) *= phi_r_symmetrized(i);

      dca::linalg::matrixop::gemm(T_r_to_k, T_k_to_r_scaled, T_k_to_k);

      for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
        for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
          T_symmetrized(i, j) = real(T_k_to_k(i, j));

      for (int i = 0; i < target_k_dmn_t::dmn_size(); i++) {
        double result = 0;

        for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
          result += T_symmetrized(i, j);

        for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
          T_symmetrized(i, j) /= result;
      }
    }
  }
}

template <typename parameters_type, typename source_k_dmn_t, typename target_k_dmn_t>
void deconvolution_routines<parameters_type, source_k_dmn_t, target_k_dmn_t>::compute_phi_inv(
    double epsilon) {
  for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
    phi_r_inv(i) = std::abs(phi_r(i)) > epsilon ? 1. / phi_r(i) : 0.;  // 1./epsilon;
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

  {
    T_k_to_k.copy_from(T_k_to_r);
    T_k_to_r_scaled.copy_from(T_k_to_r);

    for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
      for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
        T_k_to_r_scaled(i, j) *= phi_r_inv(i);

    dca::linalg::matrixop::gemm(T_r_to_k, T_k_to_r_scaled, T_k_to_k);

    for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
      for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
        T_eps(i, j) = real(T_k_to_k(i, j));
  }
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

  T_k_to_r_scaled.copy_from(T_k_to_r);

  for (int j = 0; j < target_k_dmn_t::dmn_size(); j++)
    for (int i = 0; i < target_k_dmn_t::dmn_size(); i++)
      T_k_to_r_scaled(i, j) *= phi_r_inv(i);

  dca::linalg::matrixop::gemm(T_r_to_k, T_k_to_r_scaled, T_eps);
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_DECONVOLUTION_ROUTINES_H
