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

#include "dca/phys/models/analytic_hamiltonians/hund_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/fe_as_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/twoband_Cu.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <class Lattice>
static constexpr bool has_non_density_interaction = false;

template <class BaseLattice>
static constexpr bool has_non_density_interaction<HundLattice<BaseLattice>> = true;

template <class BaseLattice>
static constexpr bool has_non_density_interaction<FeAsLattice<BaseLattice>> = true;

template <class PointGroup>
static constexpr bool has_non_density_interaction<TwoBandCu<PointGroup>> = true;

template <class Lattice, class HType, class Parameters>
std::enable_if_t<has_non_density_interaction<Lattice>> initializeNonDensityInteraction(
    HType& interaction, const Parameters& pars) {
  Lattice::initializeNonDensityInteraction(interaction, pars);
}

template <class Lattice, class HType, class Parameters>
std::enable_if_t<has_non_density_interaction<Lattice>> initializeNonDensityInteraction(
    std::unique_ptr<HType>& interaction, const Parameters& pars) {
  interaction = std::make_unique<HType>();
  Lattice::initializeNonDensityInteraction(*interaction, pars);
}

template <class Lattice, class HType, class Parameters>
std::enable_if_t<!has_non_density_interaction<Lattice>> initializeNonDensityInteraction(
    HType& /*interaction*/, const Parameters& /*pars*/) {}

}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_TRAITS_HPP
