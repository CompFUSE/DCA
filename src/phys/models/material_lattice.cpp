#include "dca/phys/models/material_hamiltonians/material_lattice.hpp"

namespace dca {
namespace phys {
namespace models {

std::string to_str(Material mat) {
  switch (mat) {
    case Material::NiO_symmetric:
      return "NiO_symmetric";
    case Material::NiO_unsymmetric:
      return "NiO_unsymmetric";
    case Material::FeSn:
      return "FeSn";
    case Material::CuO2:
      return "CuO2";
    case Material::CuO2_1band:
      return "CuO2_1band";
    case Material::SrVO3:
      return "SrVO3";
    default:
      return "";
  }
}

template class material_lattice<Material::NiO_symmetric, dca::phys::domains::no_symmetry<3>>;
template class material_lattice<Material::NiO_unsymmetric, dca::phys::domains::no_symmetry<3>>;
template class material_lattice<Material::FeSn, dca::phys::domains::no_symmetry<3>>;

}  // namespace models
}  // namespace phys
}  // namespace dca
