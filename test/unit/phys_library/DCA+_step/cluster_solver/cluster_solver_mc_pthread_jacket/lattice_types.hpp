#ifndef LATTICE_TYPES_HPP
#define LATTICE_TYPES_HPP
// INTERNAL QUESTION: how can we make sure that the copied file is included and
// not this file?
#include "phys_library/domains/cluster/cluster_typedefs.hpp"
// choiche specific
#include "phys_library/domains/cluster/symmetries/point_groups/2D/2D_square.h"
#include "phys_library/parameters/models/analytic_Hamiltonians/lattices/2D_square_lattice.h"

namespace types {
// group symmetry
using DCA_point_group_type = D4;
// using DCA_point_group_type = no_symmetry<2>;
// using DCA_point_group_type = no_symmetry<3>;
// using DCA_point_group_type = symmetry_package_spg_lib;
// using DCA_point_group_type = S4;
// using DCA_point_group_type = S4_plus;
// using DCA_point_group_type = C4;
// using DCA_point_group_type = D4_trial;
// using DCA_point_group_type = C3;
// using DCA_point_group_type = S3;
// using DCA_point_group_type = D3;
// using DCA_point_group_type = C6;
// using DCA_point_group_type = S6;
// using DCA_point_group_type = D6;
// using DCA_point_group_type = Oh;
// using DCA_point_group_type = spg_symmetry_package;

// lattice type
using lattice_type = square_lattice<DCA_point_group_type>;
// using lattice_type = triangular_lattice<DCA_point_group_type>;
// using lattice_type = bilayer_lattice<DCA_point_group_type>;
// using lattice_type = cubic_lattice<DCA_point_group_type>;
// using lattice_type = material_lattice<NiO,DCA_point_group_type>;
// using lattice_type = material_lattice<CuO2,DCA_point_group_type>;
// typedef Bett_cluster_square_2D  <DCA_point_group_type> lattice_type;
// using lattice_type = cuprate_three_band_LATTICE<DCA_point_group_type>;
// using lattice_type = Bett_cluster_square_2D_singlets<DCA_point_group_type>;
// using lattice_type = Bett_cluster_square_3D<DCA_point_group_type>;
// using lattice_type = Bett_cluster_triangular_2D<DCA_point_group_type>;
// using lattice_type = Andersen_Hamiltonians<DCA_point_group_type>;
// using lattice_type = cuprate_LATTICE_3_bands<DCA_point_group_type>;
// using lattice_type = Bett_cluster_square_2D_layered<DCA_point_group_type>;
// using lattice_type =
// Bett_cluster_square_2D_layered_non_interacting<DCA_point_group_type>;
// using lattice_type = bipartite_square_LATTICE<DCA_point_group_type>;
// using lattice_type = cuprate_2_band_LATTICE<DCA_point_group_type>;
// using lattice_type = cuprate_3_band_LATTICE<DCA_point_group_type>;
// using lattice_type = cuprate_3_band_LATTICE_kent<DCA_point_group_type>;
// using lattice_type = non_interacting_cluster_square_2D<DCA_point_group_type>;

}  // namespace types

#endif  // LATTICE_TYPES_HPP
