// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides a material lattice class that is templated on the material type.

#ifndef DCA_PHYS_MODELS_MATERIAL_HAMILTONIANS_MATERIAL_LATTICE_HPP
#define DCA_PHYS_MODELS_MATERIAL_HAMILTONIANS_MATERIAL_LATTICE_HPP

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/csv/csv_reader.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::
  
//this can be used as an element in a std::conjunction, std::disjunction
template <typename T>
using Is3D_t = typename std::conditional_t<T::DIM == 3, std::true_type, std::false_type>;

template <typename T>
using Is3D = std::enable_if_t<Is3D_t<T>::value>;

template <typename T>
using Is2D = std::enable_if_t<std::is_floating_point<T>::value, bool>;

enum class Material { NiO_symmetric, NiO_unsymmetric, CuO2, CuO2_1band, SrVO3, FeSn };
  
// forward declaration of class template
// we don't make this empty because we want this to fail at compile time if a particle specialization is unavailable.
template <Material MATERIAL, typename point_group_type, typename Enable = void>
class material_lattice;

// // Specialization for CuO2
#include "CuO2/material_lattice_CuO2.inc"

// // Specialization for CuO2 1-band
#include "CuO2_1band/material_lattice_CuO2_1band.inc"

// Specialization for NiO
#include "NiO/material_lattice_NiO.inc"

// // Specialization for SrVO3
#include "SrVO3/material_lattice_SrVO3.inc"

// // Specialization for SrVO3
#include "FeSn/material_lattice_FeSn.inc"

std::string to_str(Material mat);
  
}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_MATERIAL_HAMILTONIANS_MATERIAL_LATTICE_HPP
