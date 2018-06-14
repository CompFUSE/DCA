// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines a tetrahedron class that is templated on its dimension.

#ifndef DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_TETRAHEDRON_HPP
#define DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_TETRAHEDRON_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

#include "dca/math/geometry/tetrahedron_mesh/tetrahedron_eigenvalue_degeneracy.hpp"
#include "dca/math/util/comparison_methods.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/util/plot.hpp"

namespace dca {
namespace math {
namespace geometry {
// dca::math::geometry::

// Empty class template.
template <int dimension>
struct tetrahedron {};

// Full template specialization for a 1D tetrahedron.
#include "tetrahedron_1d.inc"

// Full template specialization for a 2D tetrahedron.
#include "tetrahedron_2d.inc"

// Full template specialization for a 3D tetrahedron.
#include "tetrahedron_3d.inc"

}  // geometry
}  // math
}  // dca

#endif  // DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_TETRAHEDRON_HPP
