// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_ROUTINES_H
#define PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_ROUTINES_H

#include "comp_library/function_library/domains/special_domains/dmn_0.h"
#include "comp_library/function_library/domains/special_domains/dmn_variadic.h"
#include "phys_library/domains/cluster/centered_cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"

namespace DCA {

template <typename parameters_type, typename source_k_dmn, typename target_k_dmn>
class interpolation_routines {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using source_r_cluster_type = typename source_k_dmn::parameter_type::dual_type;
  using r_centered_dmn = dmn_0<centered_cluster_domain<source_r_cluster_type>>;

  using w = dmn_0<frequency_domain>;

  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index

  using nu_nu_r_centered = dmn_variadic<nu, nu, r_centered_dmn>;
  using nu_nu_r_centered_w = dmn_variadic<nu, nu, r_centered_dmn, w>;

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
}

#endif  // PHYS_LIBRARY_DCA_STEP_LATTICE_MAPPING_INTERPOLATION_INTERPOLATION_ROUTINES_H
