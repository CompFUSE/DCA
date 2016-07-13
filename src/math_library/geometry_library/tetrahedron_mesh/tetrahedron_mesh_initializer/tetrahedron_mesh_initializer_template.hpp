// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Empty class declaration for the tetrahedron mesh initializer object that is templated over the
// dimension of the tetrahedron mesh and the momentum space cluster type.

#ifndef MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_MESH_INITIALIZER_TETRAHEDRON_MESH_INITIALIZER_TEMPLATE_HPP
#define MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_MESH_INITIALIZER_TETRAHEDRON_MESH_INITIALIZER_TEMPLATE_HPP

namespace math_algorithms {

template <int DIMENSION, class k_cluster_type>
class tetrahedron_mesh_initializer {};
}

#endif  // MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_MESH_INITIALIZER_TETRAHEDRON_MESH_INITIALIZER_TEMPLATE_HPP
