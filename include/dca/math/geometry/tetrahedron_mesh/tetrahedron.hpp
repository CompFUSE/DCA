// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
#include "dca/math/util/vector_operations.hpp"

#include "comp_library/function_plotting/include_plotting.h"
#include "math_library/static_functions.h"

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
