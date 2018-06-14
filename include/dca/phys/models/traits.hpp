// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE for terms of usage./
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Defines the presence of non density-density interactions using SFINAE.

#ifndef DCA_PHYS_MODELS_TRAITS_HPP
#define DCA_PHYS_MODELS_TRAITS_HPP

#include <memory>
#include <type_traits>

// INTERNAL: include with CT-INT
//#include "dca/phys/models/analytic_hamiltonians/hund_lattice.hpp"
//#include "dca/phys/models/analytic_hamiltonians/fe_as_lattice.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <class Lattice>
struct has_non_density_interaction {
  constexpr static bool value = false;
};

// INTERNAL: include with CT-INT
// template <class BaseLattice>
// struct has_non_density_interaction<HundLattice<BaseLattice>> {
//  constexpr static bool value = true;
//};
//
// template <class BaseLattice>
// struct has_non_density_interaction<FeAsLattice<BaseLattice>> {
//  constexpr static bool value = true;
//};

template <class Lattice, class HType, class Parameters>
std::enable_if_t<has_non_density_interaction<Lattice>::value, void> initializeNonDensityInteraction(
    HType& interaction, const Parameters& pars) {
  Lattice::initializeNonDensityInteraction(interaction, pars);
}

template <class Lattice, class HType, class Parameters>
std::enable_if_t<has_non_density_interaction<Lattice>::value, void> initializeNonDensityInteraction(
    std::unique_ptr<HType>& interaction, const Parameters& pars) {
  interaction.reset(new HType);
  Lattice::initializeNonDensityInteraction(*interaction, pars);
}

template <class Lattice, class HType, class Parameters>
typename std::enable_if_t<not has_non_density_interaction<Lattice>::value, void> initializeNonDensityInteraction(
    HType& /*interaction*/, const Parameters& /*pars*/) {}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_TRAITS_HPP
