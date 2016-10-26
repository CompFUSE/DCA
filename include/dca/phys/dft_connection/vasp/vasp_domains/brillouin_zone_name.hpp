// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the Brillouin zone names.

#ifndef DCA_PHYS_DFT_CONNECTION_VASP_VASP_DOMAINS_BRILLOUIN_ZONE_NAME_HPP
#define DCA_PHYS_DFT_CONNECTION_VASP_VASP_DOMAINS_BRILLOUIN_ZONE_NAME_HPP

namespace dca {
namespace phys {
namespace dft {
namespace vasp {
// dca::phys::dft::vasp::

enum brillouin_zone_name {
  NO_NAME,
  FERMI_SURFACE_SQUARE_2D_LATTICE,
  SQUARE_2D_LATTICE,
  BODY_CENTERED_TETRAGONAL_A,
  BODY_CENTERED_TETRAGONAL_B,
  SIMPLE_TETRAGONAL,
  TRICLINIC,
  FACE_CENTERED_CUBIC,
  BODY_CENTERED_CUBIC,
  SIMPLE_CUBIC,
  HEXAGONAL,
  RHOMBOHEDRAL_A,
  RHOMBOHEDRAL_B,
  SIMPLE_MONOCLINIC,
  ONE_FACE_CENTERED_MONOCLINIC_A,
  ONE_FACE_CENTERED_MONOCLINIC_B,
  SIMPLE_ORTHOROMBIC,
  BASE_CENTERED_ORTHORHOMBIC,
  BODY_CENTERED_ORTHOROMBIC,
  ALL_FACE_CENTERED_ORTHORHOMBIC_A,
  ALL_FACE_CENTERED_ORTHORHOMBIC_B
};

}  // vasp
}  // dft
}  // phys
}  // dca

#endif  // DCA_PHYS_DFT_CONNECTION_VASP_VASP_DOMAINS_BRILLOUIN_ZONE_NAME_HPP
