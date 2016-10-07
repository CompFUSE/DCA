// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Bart Ydens
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//
// This class implements the insertion and removal of segments \f$(c^{\dagger}-c\ pair)\f$.
//
// \f{eqnarray}{
//   W_{acc}(c_k \rightarrow c_{k+1}) = min \left(1, \frac{\beta l_{max}}{k+1} \frac{det F(k+1)}{det
//   F(k)}\frac{W_{loc}(k+1)}{W_{loc}(k)} \right)
// \f}

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_WALKER_TOOLS_SEGMENT_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_WALKER_TOOLS_SEGMENT_H

#include "dca/function/domains/dmn_0.hpp"

#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"

namespace QMC {

template <typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
class segment {
public:
  typedef typename configuration_t::orbital_configuration_t orbital_configuration_t;

  typedef segment<configuration_t, parameters_t, MOMS_t, concurrency_t> THIS_TYPE;

  using b = func::dmn_0<electron_band_domain>;
  using s = func::dmn_0<electron_spin_domain>;

public:
  segment(configuration_t& configuration, parameters_t& parameters, MOMS_t& MOMS,
          concurrency_t& concurrency);

  template <typename function_type_1, typename function_type_2>
  void insert_segment(int this_flavor, double mu, double& sign, function_type_1& M,
                      function_type_2& F);

  template <typename function_type_1, typename function_type_2>
  void remove_segment(int this_flavor, double mu, double& sign, function_type_1& M,
                      function_type_2& F);

private:
  configuration_t& configuration;

  concurrency_t& concurrency;
  parameters_t& parameters;
  MOMS_t& MOMS;

  int FLAVORS;
  double BETA;
};

template <typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
segment<configuration_t, parameters_t, MOMS_t, concurrency_t>::segment(
    configuration_t& configuration_ref, parameters_t& parameters_ref, MOMS_t& MOMS_ref,
    concurrency_t& concurrency_ref)
    : configuration(configuration_ref),
      concurrency(concurrency_ref),
      parameters(parameters_ref),
      MOMS(MOMS_ref) {
  FLAVORS = b::dmn_size() * s::dmn_size();
  BETA = parameters.get_beta();
}

template <typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
template <typename function_type_1, typename function_type_2>
void segment<configuration_t, parameters_t, MOMS_t, concurrency_t>::insert_segment(
    int this_flavor, double mu, double& sign, function_type_1& M, function_type_2& F) {}

template <typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
template <typename function_type_1, typename function_type_2>
void segment<configuration_t, parameters_t, MOMS_t, concurrency_t>::remove_segment(
    int this_flavor, double mu, double& sign, function_type_1& M, function_type_2& F) {}

}  // QMC

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_WALKER_TOOLS_SEGMENT_H
