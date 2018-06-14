// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the tetrahedron eigenvalue degeneracies.

#ifndef DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_TETRAHEDRON_EIGENVALUE_DEGENERACY_HPP
#define DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_TETRAHEDRON_EIGENVALUE_DEGENERACY_HPP

namespace dca {
namespace math {
namespace geometry {
// dca::math::geometry::

enum TetrahedronEigenvalueDegeneracy {
  NO_DEGENERACY,
  TWOFOLD_DEGENERACY,
  THREEFOLD_DEGENERACY,
  THREEFOLD_DEGENERACY_A,
  THREEFOLD_DEGENERACY_B,
  FOURFOLD_DEGENERACY
};

}  // geometry
}  // math
}  // dca

#endif  // DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_TETRAHEDRON_EIGENVALUE_DEGENERACY_HPP
