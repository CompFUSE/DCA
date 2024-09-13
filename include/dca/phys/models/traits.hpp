// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
// See LICENSE for terms of usage./
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// Defines the presence of non density-density interactions using SFINAE.

#ifndef DCA_PHYS_MODELS_TRAITS_HPP
#define DCA_PHYS_MODELS_TRAITS_HPP

#include <memory>
#include <type_traits>

#include "dca/function/function.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

using NuDmn = func::dmn_variadic<func::dmn_0<domains::electron_band_domain>,
                                 func::dmn_0<domains::electron_spin_domain>>;

template <typename Scalar, class Parameters>
using NonDensityIntHamiltonian =
    func::function<Scalar,
                   func::dmn_variadic<NuDmn, NuDmn, NuDmn, NuDmn, typename Parameters::RClusterDmn>>;

// // Class to detect if class T implements the templated "initializeNonDensityInteraction" method.
// template <class PARAMS>
// class HasInitializeNonDensityInteractionMethod_t {
// private:
//   template <typename PARS>
//   static auto test(const PARS& pars) ->
//   decltype(PARS::lattice_type::initializeNonDensityInteraction(std::declval<NonDensityIntHamiltonian<double,
//   PARS>>(), pars), std::true_type{});

//   //  static std::false_type test(...);

// public:
//   using type = std::is_same<decltype(test(std::declval<PARAMS>())), std::false_type>;
//   //decltype(test(std::declval<PARAMS>()));
// };

// Class to detect if class T implements the templated "initializeNonDensityInteraction" method.
template <class PARAMS>
class HasInitializeNonDensityInteractionMethod {
private:
  template <typename PARS>
  static auto test(const PARS& pars) -> decltype(PARS::lattice_type::initializeNonDensityInteraction(
      std::declval<NonDensityIntHamiltonian<typename PARS::Scalar, PARS>& >(), pars));  //, std::true_type{});

  static std::false_type test(...);

public:
  constexpr static bool value =
      !std::is_same_v<decltype(test(std::declval<PARAMS>())), std::false_type>;
};

template <class Parameters>
std::enable_if_t<HasInitializeNonDensityInteractionMethod<Parameters>::value> initializeNonDensityInteraction(
    NonDensityIntHamiltonian<typename Parameters::Scalar, Parameters>& interaction,
    const Parameters& pars) {
  Parameters::lattice_type::initializeNonDensityInteraction(interaction, pars);
}

template <typename Scalar, class Lattice, class Parameters>
std::enable_if_t<HasInitializeNonDensityInteractionMethod<Parameters>::value> initializeNonDensityInteraction(
    std::unique_ptr<NonDensityIntHamiltonian<Scalar, Parameters>>& interaction,
    const Parameters& pars) {
  interaction = std::make_unique<NonDensityIntHamiltonian<Scalar, Parameters>>();
  Parameters::lattice_type::initializeNonDensityInteraction(*interaction, pars);
}

template <class Lattice, class HType, class Parameters>
std::enable_if_t<!HasInitializeNonDensityInteractionMethod<Parameters>::value> initializeNonDensityInteraction(
    HType& /*interaction*/, const Parameters& /*pars*/) {}

}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_TRAITS_HPP
