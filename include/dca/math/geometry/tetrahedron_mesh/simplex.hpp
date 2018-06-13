// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class represents a simplex (= edge-corner) of the Brillouin zone. It is templated on the
// dimension of the k-space.

#ifndef DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_SIMPLEX_HPP
#define DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_SIMPLEX_HPP

#include <vector>

namespace dca {
namespace math {
namespace geometry {
// dca::math::geometry::

template <int dimension>
struct simplex {
  std::vector<double> k_vec;
};

}  // geometry
}  // math
}  // dca

#endif  // DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_SIMPLEX_HPP
