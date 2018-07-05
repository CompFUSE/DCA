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
// This file provides functions to compute the Gaussian quadrature weights and abscissas for 1D, 2D
// and 3D. It makes use of the SIMPLEX_GM_RULE C++ library
// (https://people.sc.fsu.edu/~jburkardt/cpp_src/simplex_gm_rule/simplex_gm_rule.html).

#ifndef DCA_MATH_GEOMETRY_GAUSSIAN_QUADRATURE_COMPUTE_WEIGHTS_AND_ABSCISSAS_HPP
#define DCA_MATH_GEOMETRY_GAUSSIAN_QUADRATURE_COMPUTE_WEIGHTS_AND_ABSCISSAS_HPP

#include "simplex_gm_rule.hpp"
#include "dca/math/geometry/tetrahedron_mesh/tetrahedron.hpp"

namespace dca {
namespace math {
namespace geometry {
// dca::math::geometry::

// Empty template declaration
template <int dimension>
void computeWeightsAndAbscissas(const int rule, tetrahedron<dimension>& tet);

// Template specializations
// 1D
template <>
void computeWeightsAndAbscissas(const int rule, tetrahedron<1>& tet);
// 2D
template <>
void computeWeightsAndAbscissas(const int rule, tetrahedron<2>& tet);
// 3D
template <>
void computeWeightsAndAbscissas(const int rule, tetrahedron<3>& tet);

}  // geometry
}  // math
}  // dca

#endif  // DCA_MATH_GEOMETRY_GAUSSIAN_QUADRATURE_COMPUTE_WEIGHTS_AND_ABSCISSAS_HPP
