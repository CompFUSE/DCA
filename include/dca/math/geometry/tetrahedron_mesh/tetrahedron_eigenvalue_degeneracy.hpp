// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
