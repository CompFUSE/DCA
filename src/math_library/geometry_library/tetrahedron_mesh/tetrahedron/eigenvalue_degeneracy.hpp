// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// This header file defines the eigenvalue degeneracies.

#ifndef MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_EIGENVALUE_DEGENERACY_HPP
#define MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_EIGENVALUE_DEGENERACY_HPP

namespace math_algorithms {

enum eigenvalue_degeneracy {
  NO_DEGENERACY,
  TWOFOLD_DEGENERACY,
  THREEFOLD_DEGENERACY,
  THREEFOLD_DEGENERACY_A,
  THREEFOLD_DEGENERACY_B,
  FOURFOLD_DEGENERACY
};

}  // math_algorithms

#endif  // MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_EIGENVALUE_DEGENERACY_HPP
