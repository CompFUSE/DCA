#ifndef MODEL_TYPES_HPP
#define MODEL_TYPES_HPP
#include "phys_library/parameters/models/tight_binding_model.h"

namespace types {

using model = tight_binding_model<lattice_type>;
// using LATTICE = dft_LATTICE<3,DCA_point_group_type>;
// using LATTICE = Koshevnikov_LATTICE;
// using LATTICE = cuprate_single_band_LATTICE<DCA_point_group_type>;
// using LATTICE = cuprate_three_band_LATTICE<DCA_point_group_type>;
// using LATTICE = bilayer_LATTICE<DCA_point_group_type>;
} // end namespace types

#endif