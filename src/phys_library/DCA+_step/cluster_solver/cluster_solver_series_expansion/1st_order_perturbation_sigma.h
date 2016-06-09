// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// This class computes the self-energy in first order.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_1ST_ORDER_PERTURBATION_SIGMA_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_1ST_ORDER_PERTURBATION_SIGMA_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/template_perturbation_sigma.h"

#include <complex>

#include "comp_library/function_library/include_function_library.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/compute_interaction.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"

namespace DCA {
namespace SERIES_EXPANSION {

template <class parameter_type, class k_dmn_t>
class sigma_perturbation<1, parameter_type, k_dmn_t> {
public:
  using U_type = typename compute_interaction::function_type;

  using w = dmn_0<frequency_domain>;
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index

  typedef FUNC_LIB::function<std::complex<double>, dmn_variadic<nu, nu, k_dmn_t, w>> function_type;

public:
  sigma_perturbation(parameter_type& parameters_ref, compute_interaction& interaction_obj);

  function_type& get_function() {
    return Sigma;
  }

  void execute_on_cluster(FUNC_LIB::function<double, nu>& occupancy);

  template <typename Writer>
  void write(Writer& writer);

protected:
  parameter_type& parameters;

  U_type& U;

  FUNC_LIB::function<std::complex<double>, dmn_variadic<nu, nu>> d_matrix;
  FUNC_LIB::function<std::complex<double>, dmn_variadic<nu, nu, k_dmn_t, w>> Sigma;
};

template <class parameter_type, class k_dmn_t>
sigma_perturbation<1, parameter_type, k_dmn_t>::sigma_perturbation(parameter_type& parameters_ref,
                                                                   compute_interaction& interaction_obj)
    : parameters(parameters_ref),

      U(interaction_obj.get_function()),

      d_matrix("d-matrix"),
      Sigma("Sigma-1st-order") {}

template <class parameter_type, class k_dmn_t>
template <typename Writer>
void sigma_perturbation<1, parameter_type, k_dmn_t>::write(Writer& /*writer*/) {}

template <class parameter_type, class k_dmn_t>
void sigma_perturbation<1, parameter_type, k_dmn_t>::execute_on_cluster(
    FUNC_LIB::function<double, nu>& occupancy) {
  Sigma = 0.;

  for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind) {
    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
      for (int i = 0; i < b::dmn_size(); ++i)
        for (int si = 0; si < s::dmn_size(); ++si)
          for (int j = 0; j < b::dmn_size(); ++j)
            for (int sj = 0; sj < s::dmn_size(); ++sj)
              Sigma(i, si, i, si, k_ind, w_ind) += U(i, si, j, sj) * (occupancy(j, sj) - 1. / 2.);
    }
  }
}
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_1ST_ORDER_PERTURBATION_SIGMA_H
