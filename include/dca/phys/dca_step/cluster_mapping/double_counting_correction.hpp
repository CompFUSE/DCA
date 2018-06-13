// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//
// This class perfoms the double counting correction.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_DOUBLE_COUNTING_CORRECTION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_DOUBLE_COUNTING_CORRECTION_HPP

#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename parameters_type, typename MOMS_type>
class double_counting_correction {
public:
  using concurrency_type = typename parameters_type::concurrency_type;

  using w = func::dmn_0<domains::frequency_domain>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using DCA_k_cluster_type =
      domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_DCA = func::dmn_0<DCA_k_cluster_type>;

public:
  double_counting_correction(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  void initialize();

  void execute_before_solver();

  void execute_after_solver();

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  MOMS_type& MOMS;

  func::function<double, nu> mu_HALF;
  func::function<double, nu> DC;
};

template <typename parameters_type, typename MOMS_type>
double_counting_correction<parameters_type, MOMS_type>::double_counting_correction(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      MOMS(MOMS_ref),

      mu_HALF("mu_HALF"),
      DC("DC") {
  initialize();
}

template <typename parameters_type, typename MOMS_type>
void double_counting_correction<parameters_type, MOMS_type>::initialize() {
  const std::vector<int>& interacting_bands = parameters.get_interacting_orbitals();

  for (int b_i_int = 0; b_i_int < interacting_bands.size(); b_i_int++) {
    int b_i = interacting_bands[b_i_int];

    for (int s_i = 0; s_i < s::dmn_size(); s_i++) {
      for (int b_j_int = 0; b_j_int < interacting_bands.size(); b_j_int++) {
        int b_j = interacting_bands[b_j_int];

        for (int s_j = 0; s_j < s::dmn_size(); s_j++)
          mu_HALF(b_i, s_i) += (1. / 2.) * (MOMS.H_interactions(b_j, s_j, b_i, s_i, 0) +
                                            MOMS.H_interactions(b_i, s_i, b_j, s_j, 0)) *
                               1. / 2.;
      }
    }
  }

  if (parameters.get_double_counting_method() == "constant-correction-without-U-correction")
    for (int b_ind = 0; b_ind < interacting_bands.size(); b_ind++)
      for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++)
        DC(interacting_bands[b_ind], s_ind) =
            parameters.get_double_counting_correction();  //-mu_HALF(interacting_bands[b_ind],s_ind);

  if (parameters.get_double_counting_method() == "constant-correction-with-U-correction")
    for (int b_ind = 0; b_ind < interacting_bands.size(); b_ind++)
      for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++)
        DC(interacting_bands[b_ind], s_ind) =
            parameters.get_double_counting_correction() - mu_HALF(interacting_bands[b_ind], s_ind);
}

template <typename parameters_type, typename MOMS_type>
void double_counting_correction<parameters_type, MOMS_type>::execute_before_solver() {
  if (parameters.get_double_counting_method() != "none") {
    const std::vector<int>& interacting_bands = parameters.get_interacting_orbitals();

    if (parameters.get_double_counting_method() == "constant-correction-with-U-correction" or
        parameters.get_double_counting_method() == "constant-correction-without-U-correction") {
      for (int b_ind = 0; b_ind < interacting_bands.size(); b_ind++)
        for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++)
          for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
            for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++)
              MOMS.Sigma_cluster(interacting_bands[b_ind], s_ind, interacting_bands[b_ind], s_ind,
                                 k_ind, w_ind) += DC(interacting_bands[b_ind], s_ind);
    }
  }
}

template <typename parameters_type, typename MOMS_type>
void double_counting_correction<parameters_type, MOMS_type>::execute_after_solver() {
  if (parameters.get_double_counting_method() != "none") {
    const std::vector<int>& interacting_bands = parameters.get_interacting_orbitals();

    if (parameters.get_double_counting_method() == "constant-correction-with-U-correction" or
        parameters.get_double_counting_method() == "constant-correction-without-U-correction") {
      for (int b_ind = 0; b_ind < interacting_bands.size(); b_ind++)
        for (int s_ind = 0; s_ind < s::dmn_size(); s_ind++)
          for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
            for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++)
              MOMS.Sigma(interacting_bands[b_ind], s_ind, interacting_bands[b_ind], s_ind, k_ind,
                         w_ind) -= DC(interacting_bands[b_ind], s_ind);
    }
  }
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_DOUBLE_COUNTING_CORRECTION_HPP
