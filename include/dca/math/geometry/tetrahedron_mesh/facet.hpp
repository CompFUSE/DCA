// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class represents a facet (= edge-surface) of the Brillouin zone. It is templated on the
// dimension of the k-space.

#ifndef DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_FACET_HPP
#define DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_FACET_HPP

#include <vector>
#include "dca/math/geometry/tetrahedron_mesh/simplex.hpp"

namespace dca {
namespace math {
namespace geometry {
// dca::math::geometry::

// Empty template class declaration
template <int DIMENSION>
struct facet {};

// Full template specialiaztion for 2D
template <>
struct facet<2> {
  std::vector<int> index;

  static void find_linear_parameters(int* coor, double* parameters,
                                     std::vector<simplex<2>>& simplex_vector);
  static bool is_facet(int* coor, std::vector<simplex<2>>& simplex_vector);
  static bool is_facet(int* coor, std::vector<std::vector<double>>& vectors);
  static bool equal(facet& f1, facet& f2, std::vector<simplex<2>>& simplex_vector);
};

// Full template specialiaztion for 3D
template <>
struct facet<3> {
  std::vector<int> index;

  static void find_linear_parameters(int* coor, double* parameters,
                                     std::vector<simplex<3>>& simplex_vector);
  static bool is_facet(int* coor, std::vector<simplex<3>>& simplex_vector);
  static bool is_facet(int* coor, std::vector<std::vector<double>>& simplex_vector);
  static bool equal(facet& f1, facet& f2, std::vector<simplex<3>>& simplex_vector);
};

}  // geometry
}  // math
}  // dca

#endif  // DCA_MATH_GEOMETRY_TETRAHEDRON_MESH_FACET_HPP
