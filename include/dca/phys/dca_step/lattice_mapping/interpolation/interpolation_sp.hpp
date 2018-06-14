// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes the interpolated cluster self-energy using the alpha transformation.

#ifndef DCA_PHYS_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_SP_HPP
#define DCA_PHYS_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_SP_HPP

#include <complex>
#include <iostream>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/phys/dca_step/lattice_mapping/interpolation/interpolation_routines.hpp"
#include "dca/phys/dca_step/lattice_mapping/interpolation/transform_to_alpha.hpp"
#include "dca/phys/domains/cluster/centered_cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace latticemapping {
// dca::phys::latticemapping::

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
class interpolation_sp : public interpolation_routines<parameters_type, source_k_dmn, target_k_dmn> {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using source_r_cluster_type = typename source_k_dmn::parameter_type::dual_type;
  using r_centered_dmn = func::dmn_0<domains::centered_cluster_domain<source_r_cluster_type>>;

  using w = func::dmn_0<domains::frequency_domain>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using DCA_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_DCA = func::dmn_0<DCA_k_cluster_type>;
  using host_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::LATTICE_SP,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_HOST = func::dmn_0<host_k_cluster_type>;

  using nu_nu_r_centered = func::dmn_variadic<nu, nu, r_centered_dmn>;
  using nu_nu_r_centered_w = func::dmn_variadic<nu, nu, r_centered_dmn, w>;
  using nu_nu_k_DCA_w = func::dmn_variadic<nu, nu, k_DCA, w>;
  using nu_nu_k_HOST_w = func::dmn_variadic<nu, nu, k_HOST, w>;

public:
  interpolation_sp(parameters_type& parameters_ref);

  void execute_with_alpha_transformation(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn, w>>& cluster_self_energy,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn, w>>& interp_self_energy);

private:
  void execute(
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn>>& cluster_self_energy,
      func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn>>& interp_self_energy);

  void execute(func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn, w>>&
                   cluster_self_energy,
               func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn, w>>&
                   interp_self_energy);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;
};

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
interpolation_sp<parameters_type, source_k_dmn, target_k_dmn>::interpolation_sp(
    parameters_type& parameters_ref)
    : interpolation_routines<parameters_type, source_k_dmn, target_k_dmn>(parameters_ref),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t" << __FUNCTION__ << " is created " << dca::util::print_time();
}

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
void interpolation_sp<parameters_type, source_k_dmn, target_k_dmn>::execute_with_alpha_transformation(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn, w>>& cluster_self_energy,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn, w>>& interp_self_energy) {
  r_centered_dmn::parameter_type::initialize();

  func::function<std::complex<double>, nu_nu_k_DCA_w> cluster_alpha_k("cluster_alpha_k");
  func::function<std::complex<double>, nu_nu_k_HOST_w> interpolated_alpha_k("interpolated_alpha");

  transform_to_alpha::forward(1., cluster_self_energy, cluster_alpha_k);

  execute(cluster_alpha_k, interpolated_alpha_k);

  transform_to_alpha::backward(1., interp_self_energy, interpolated_alpha_k);
}

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
void interpolation_sp<parameters_type, source_k_dmn, target_k_dmn>::execute(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn, w>>& cluster_function,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn, w>>& interp_function) {
  r_centered_dmn::parameter_type::initialize();

  func::function<std::complex<double>, nu_nu_r_centered_w> cluster_centered_function(
      "cluster_centered_function");

  math::transform::FunctionTransform<source_k_dmn, r_centered_dmn>::execute(
      cluster_function, cluster_centered_function);

  math::transform::FunctionTransform<r_centered_dmn, target_k_dmn>::execute(
      cluster_centered_function, interp_function);
}

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
void interpolation_sp<parameters_type, source_k_dmn, target_k_dmn>::execute(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, source_k_dmn>>& cluster_function,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, target_k_dmn>>& interp_function) {
  r_centered_dmn::parameter_type::initialize();

  func::function<std::complex<double>, nu_nu_r_centered> cluster_centered_function(
      "cluster_centered_function");

  math::transform::FunctionTransform<source_k_dmn, r_centered_dmn>::execute(
      cluster_function, cluster_centered_function);

  math::transform::FunctionTransform<r_centered_dmn, target_k_dmn>::execute(
      cluster_centered_function, interp_function);
}

}  // latticemapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_SP_HPP
