//-*-C++-*-

#ifndef MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_SIMPLEX_H
#define MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_SIMPLEX_H

namespace math_algorithms {
/*!
 *  \class   simplex
 *  \ingroup TETRAHEDRON
 *
 *  \author Peter Staar
 *  \brief  This class represents a simplex (= edge-corner) of the Brillouin-zone. It is templated
 * over the dimension of k-space.
 */
template <int dimension>
struct simplex {
public:
  std::vector<double> k_vec;
};
}

#endif  // MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_SIMPLEX_H
