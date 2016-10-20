// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef PHYS_LIBRARY_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_BRILLOUIN_ZONE_CUT_DOMAIN_HPP
#define PHYS_LIBRARY_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_BRILLOUIN_ZONE_CUT_DOMAIN_HPP

#include <cassert>
#include <string>
#include <vector>

#include "dca/phys/domains/cluster/cluster_specifications.hpp"

namespace DFT {
namespace VASP {
// DFT::VASP::

enum brillouin_zone_names {
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

typedef brillouin_zone_names brillouin_zone_name;

class vasp_brillouin_zone_cut_domain {
  const static int DIMENSION = 3;

public:
  typedef domains::cluster_specifications<double, domains::VASP_LATTICE, domains::MOMENTUM_SPACE,
                                          domains::PARALLELLEPIPEDUM>::dmn_specifications_type
      dmn_specifications_type;

  typedef vasp_brillouin_zone_cut_domain this_type;
  typedef std::vector<double> element_type;

public:
  static int& get_N();

  static brillouin_zone_name& get_brillouin_zone_name();

  static std::string& get_name();

  static int& get_size();
  static std::vector<element_type>& get_elements();

  static void initialize(std::vector<std::vector<double>>& path_vecs);
};

int& vasp_brillouin_zone_cut_domain::get_N() {
  static int N = 65;
  return N;
}

brillouin_zone_name& vasp_brillouin_zone_cut_domain::get_brillouin_zone_name() {
  static brillouin_zone_name name = NO_NAME;
  return name;
}

std::string& vasp_brillouin_zone_cut_domain::get_name() {
  static std::string name = 0;
  return name;
}

int& vasp_brillouin_zone_cut_domain::get_size() {
  static int size = 0;
  return size;
}

std::vector<std::vector<double>>& vasp_brillouin_zone_cut_domain::get_elements() {
  static std::vector<std::vector<double>> v(0);
  return v;
}

void vasp_brillouin_zone_cut_domain::initialize(std::vector<std::vector<double>>& path_vecs) {
  get_elements().resize(0);

  std::vector<double> k(DIMENSION);

  for (int i = 0; i < int(path_vecs.size()) - 1; i++) {
    std::vector<double>& p0 = path_vecs[i];
    std::vector<double>& p1 = path_vecs[i + 1];

    assert(DIMENSION == int(p0.size()));
    assert(DIMENSION == int(p1.size()));

    for (int l = 0; l < get_N(); l++) {
      for (int z = 0; z < DIMENSION; z++)
        k[z] = (1. - double(l) / double(get_N())) * p0[z] + double(l) / double(get_N()) * p1[z];

      get_elements().push_back(k);
    }
  }

  get_size() = get_elements().size();
}

}  // VASP
}  // DFT

#endif  // PHYS_LIBRARY_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_BRILLOUIN_ZONE_CUT_DOMAIN_HPP
