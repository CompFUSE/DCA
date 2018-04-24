// Copyright (C) 2018 ETH Zurich
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

#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/twoband_chain.hpp"
#include "dca/phys/models/analytic_hamiltonians/singleband_chain.hpp"
// #include "dca/phys/models/analytic_hamiltonians/fourband_lattice.hpp"
// #include "dca/phys/models/analytic_hamiltonians/twoband_lattice.hpp"
#include "dca/phys/models/material_hamiltonians/material_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

// Empty class template
template <typename Model>
class ModelParameters {};

// Specialization for 2D 2-band model
// #include "model_parameters_2d_2band.inc"

// Specialization for 2D 4-band model
// #include "model_parameters_2d_4band.inc"

// Specialization for square lattice bilayer Hubbard model
#include "model_parameters_bilayer_hubbard.inc"

// Specialization for material models
#include "model_parameters_material.inc"

// Specialization for single-band Hubbard model
#include "model_parameters_single_band_hubbard.inc"

#include "model_parameters_singleband_chain.inc"
#include "model_parameters_twoband_chain.inc"

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_MODEL_PARAMETERS_HPP
