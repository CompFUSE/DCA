#include "dca/phys/models/material_hamiltonians/material_lattice.hpp"

namespace dca {
namespace phys {
namespace models {
template class material_lattice<Material::NiO_symmetric, dca::phys::domains::no_symmetry<3>>;
template class material_lattice<Material::NiO_unsymmetric, dca::phys::domains::no_symmetry<3>>;
template class material_lattice<Material::FeSn, dca::phys::domains::no_symmetry<3>>;
}  // namespace models
}  // namespace phys
}  // namespace dca
