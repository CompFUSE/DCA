// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides utility functions for tasks that are common between different lattices.

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_UTIL_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_UTIL_HPP

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"

namespace dca {
namespace phys {
namespace models {
namespace util {
// dca::phys::models::detail::

// Initializes the interaction part of the single-band real space Hubbard Hamiltonian with on-site
// and nearest-neighbor interaction.
// nn_vec contains the vectors that define the nearest neighbors.
// In: parameters
//     nn_vec
// Out: H_int
template <typename ParametersType, typename BandDmn, typename SpinDmn, typename RDmn>
void initializeSingleBandHint(
    const ParametersType& parameters,
    const std::vector<typename RDmn::parameter_type::element_type>& nn_vec,
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_int) {
  if (BandDmn::dmn_size() != 1)
    throw std::logic_error("Band domain size must be 1.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  // Get the index of the origin (0,0).
  const int origin = RDmn::parameter_type::origin_index();

  // Compute indices of nearest neighbors (nn) w.r.t. origin.
  std::vector<int> nn_index;

  const std::vector<typename RDmn::parameter_type::element_type>& super_basis =
      RDmn::parameter_type::get_super_basis_vectors();
  const std::vector<typename RDmn::parameter_type::element_type>& elements =
      RDmn::parameter_type::get_elements();

  for (const auto& vec : nn_vec) {
    std::vector<double> nn_vec_translated =
        domains::cluster_operations::translate_inside_cluster(vec, super_basis);
    nn_index.push_back(
        domains::cluster_operations::index(nn_vec_translated, elements, domains::BRILLOUIN_ZONE));
  }

  // Set all elements to zero.
  H_int = 0.;

  // Nearest-neighbor opposite spin interaction
  const double V = parameters.get_V();
  for (auto index : nn_index) {
    H_int(0, 0, 0, 1, index) = V;
    H_int(0, 1, 0, 0, index) = V;
  }

  // Nearest-neighbor same spin interaction
  const double V_prime = parameters.get_V_prime();
  for (auto index : nn_index) {
    H_int(0, 0, 0, 0, index) = V_prime;
    H_int(0, 1, 0, 1, index) = V_prime;
  }

  // On-site interaction
  // This has to be set last since for small clusters a nearest neighbor might
  // be the same site and therefore V would overwrite U.
  const double U = parameters.get_U();
  H_int(0, 0, 0, 1, origin) = U;
  H_int(0, 1, 0, 0, origin) = U;
}

}  // util
}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_UTIL_HPP
