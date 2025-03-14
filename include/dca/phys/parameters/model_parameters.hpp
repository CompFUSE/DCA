// Copyright (C) 2018 ETH ZurichOB
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class reads, stores, and writes the model paramters.
// It is templated on the model type and only implemented when specialized for the bilayer Hubbard
// model, single-band Hubbard model (templated on the type of the lattice) and material models.

#ifndef DCA_PHYS_PARAMETERS_MODEL_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_MODEL_PARAMETERS_HPP

#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/La3Ni2O7_bilayer.hpp"
#include "dca/phys/models/analytic_hamiltonians/twoOrbital.hpp"
#include "dca/phys/models/analytic_hamiltonians/fe_as_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/twoband_chain.hpp"
#include "dca/phys/models/analytic_hamiltonians/singleband_chain.hpp"
#include "dca/phys/models/analytic_hamiltonians/threeband_hubbard.hpp"
// #include "dca/phys/models/analytic_hamiltonians/fourband_lattice.hpp"
// #include "dca/phys/models/analytic_hamiltonians/twoband_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/Kagome_hubbard.hpp"
#include "dca/phys/models/material_hamiltonians/material_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/hund_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/rashba_hubbard.hpp"
#include "dca/phys/models/analytic_hamiltonians/Moire_Hubbard.hpp"
#include "dca/phys/models/analytic_hamiltonians/twoband_Cu.hpp"
#include "dca/phys/models/analytic_hamiltonians/Kagome_hubbard.hpp"
#include "dca/phys/models/tight_binding_model.hpp"

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

// Empty class template
template <typename Model>
class ModelParameters {};

template <typename T, typename = void>
struct HasCustomSpin : std::false_type {};

template <typename T>
struct HasCustomSpin<T, decltype(std::declval<T>().SPIN, void())> : std::true_type{};

// Specialization for 2D 2-band model
// #include "model_parameters_2d_2band.inc"

// Specialization for 2D 4-band model
// #include "model_parameters_2d_4band.inc"

#include "model_parameters_Kagome_hubbard.inc"

// Specialization for square lattice bilayer Hubbard model
#include "model_parameters_bilayer_hubbard.inc"

// Specialization for square lattice two-orbital bilayer Hubbard model for La3Ni2O7
#include "model_parameters_La3Ni2O7_bilayer.inc"

// Specialization for generic square lattice two-orbital model
#include "model_parameters_twoOrbital.inc"

// Specialization for FeAs superconducting model.
#include "model_parameters_fe_as.inc"

// Specialization for 2D bilayer model with a spin flip term
#include "model_parameters_hund.inc"

// Specialization for material models
#include "model_parameters_material.inc"

// Specialization for single-band Hubbard model
#include "model_parameters_single_band_hubbard.inc"

// Specialization for twoband Cu model
#include "twoband_Cu_parameters.inc"

// Specialization for Rashba-Hubbard model
#include "model_parameters_rashba_hubbard.inc"

// Specialization for Moire-Hubbard model
#include "model_parameters_moire_hubbard.inc"

#include "model_parameters_singleband_chain.inc"
#include "model_parameters_twoband_chain.inc"
#include "model_parameters_threeband_hubbard.inc"

}  // namespace params
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_PARAMETERS_MODEL_PARAMETERS_HPP
