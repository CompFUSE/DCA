#ifndef LATTICE_TYPES_HPP
#define LATTICE_TYPES_HPP
// INTERNAL QUESTION: how can we make sure that the copied file is included and not this file?
#include "phys_library/domains/cluster/cluster_typedefs.hpp"
// choiche specific
#include "phys_library/parameters/models/analytic_Hamiltonians/lattices/2D_bilayer_lattice.h"
#include "phys_library/domains/cluster/symmetries/point_groups/2D/2D_square.h"

namespace types {
// group symmetry
    typedef D4               DCA_point_group_type;

// lattice type
    typedef bilayer_lattice<DCA_point_group_type> lattice_type;

}  // namespace types

#endif  // LATTICE_TYPES_HPP
