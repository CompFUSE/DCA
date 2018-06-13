// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class provides routines for the interpolation step.
//
// TODO: This class doesn't do anything except for providing type definitions.

#ifndef DCA_PHYS_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_ROUTINES_HPP
#define DCA_PHYS_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_ROUTINES_HPP

#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/phys/domains/cluster/centered_cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

namespace dca {
namespace phys {
namespace latticemapping {
// dca::phys::latticemapping::

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
class interpolation_routines {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using source_r_cluster_type = typename source_k_dmn::parameter_type::dual_type;
  using r_centered_dmn = func::dmn_0<domains::centered_cluster_domain<source_r_cluster_type>>;

  using w = func::dmn_0<domains::frequency_domain>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using nu_nu_r_centered = func::dmn_variadic<nu, nu, r_centered_dmn>;
  using nu_nu_r_centered_w = func::dmn_variadic<nu, nu, r_centered_dmn, w>;

public:
  interpolation_routines(parameters_type& parameters_ref);

private:
  void initialize();

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  // MATH_LIBRARY::gaussian_fit<double, source_k_dmn, target_k_dmn> gaussian_fit_obj;
};

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
interpolation_routines<parameters_type, source_k_dmn, target_k_dmn>::interpolation_routines(
    parameters_type& parameters_ref)
    : parameters(parameters_ref), concurrency(parameters.get_concurrency()) {
  initialize();
}

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
void interpolation_routines<parameters_type, source_k_dmn, target_k_dmn>::initialize() {
  // gaussian_fit_obj.initialize_K_to_k(true, 1.e-3);
}

}  // latticemapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_ROUTINES_HPP
