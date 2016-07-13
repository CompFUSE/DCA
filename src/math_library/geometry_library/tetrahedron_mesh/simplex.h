// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class represents a simplex (= edge-corner) of the Brillouin-zone. It is templated over the
// dimension of k-space.

#ifndef MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_SIMPLEX_H
#define MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_SIMPLEX_H

#include <vector>

namespace math_algorithms {

template <int dimension>
struct simplex {
public:
  std::vector<double> k_vec;
};

}  // math_algorithms

#endif  // MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_SIMPLEX_H
