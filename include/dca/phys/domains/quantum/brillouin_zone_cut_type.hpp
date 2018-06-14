// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file defines the various cuts of the Brillouin zone.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_BRILLOUIN_ZONE_CUT_TYPE_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_BRILLOUIN_ZONE_CUT_TYPE_HPP

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

enum BrillouinZoneCutType {
  BRILLOUIN_ZONE_CUT_TEMPLATE,
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

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_BRILLOUIN_ZONE_CUT_TYPE_HPP
