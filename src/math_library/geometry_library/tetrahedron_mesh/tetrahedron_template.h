//-*-C++-*-

#ifndef MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_TEMPLATE_H
#define MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_TEMPLATE_H
enum eigenvalue_degeneracy {
  NO_DEGENERACY,
  TWOFOLD_DEGENERACY,
  THREEFOLD_DEGENERACY,
  THREEFOLD_DEGENERACY_A,
  THREEFOLD_DEGENERACY_B,
  FOURFOLD_DEGENERACY
};

namespace math_algorithms {
/*! \defgroup TETRAHEDRON
 *  \ingroup  TETRAHEDRON-MESH
 */

/*! \class   tetrahedron
 *  \ingroup TETRAHEDRON
 *
 *  \author  Peter Staar
 *  \brief   Empty template for a tetrahedron object, templated over its dimension.
 */

template <int dimension>
struct tetrahedron {};
}

#endif  // MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_H
