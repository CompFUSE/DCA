// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements brillouin_zone_path_domain.hpp.

#include "dca/phys/domains/quantum/brillouin_zone_path_domain.hpp"
#include <cmath>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

std::vector<std::vector<double>> brillouin_zone_path_domain<SQUARE_2D_LATTICE>::initialize_elements() {
  std::vector<element_type> collection_k_vecs(0);

  std::vector<double> b0(2, 0);
  std::vector<double> b1(2, 0);

  b0[0] = 2 * M_PI;
  b0[1] = 0;
  b1[0] = 0;
  b1[1] = 2 * M_PI;

  std::vector<double> k0(2);
  std::vector<double> k1(2);
  std::vector<double> k2(2);
  std::vector<double> k3(2);

  for (int i = 0; i < 2; i++) {
    k0[i] = 0 * b0[i] + 0. * b1[i];
    k1[i] = 1. / 2. * b0[i] + 1. / 2. * b1[i];
    k2[i] = 1. / 2. * b0[i] + 0. * b1[i];
    k3[i] = 0. * b0[i] + 0. * b1[i];
  }

  for (int l = 0; l < INTERPOLATION_NB; l++) {
    std::vector<double> k(2, 0);

    k[0] = (1. - double(l) / double(INTERPOLATION_NB)) * k0[0] +
           double(l) / double(INTERPOLATION_NB) * k1[0];
    k[1] = (1. - double(l) / double(INTERPOLATION_NB)) * k0[1] +
           double(l) / double(INTERPOLATION_NB) * k1[1];

    collection_k_vecs.push_back(k);
  }

  for (int l = 0; l < INTERPOLATION_NB; l++) {
    std::vector<double> k(2, 0);

    k[0] = (1. - double(l) / double(INTERPOLATION_NB)) * k1[0] +
           double(l) / double(INTERPOLATION_NB) * k2[0];
    k[1] = (1. - double(l) / double(INTERPOLATION_NB)) * k1[1] +
           double(l) / double(INTERPOLATION_NB) * k2[1];

    collection_k_vecs.push_back(k);
  }

  for (int l = 0; l < INTERPOLATION_NB; l++) {
    std::vector<double> k(2, 0);

    k[0] = (1. - double(l) / double(INTERPOLATION_NB)) * k2[0] +
           double(l) / double(INTERPOLATION_NB) * k3[0];
    k[1] = (1. - double(l) / double(INTERPOLATION_NB)) * k2[1] +
           double(l) / double(INTERPOLATION_NB) * k3[1];

    collection_k_vecs.push_back(k);
  }

  collection_k_vecs.push_back(k3);

  return collection_k_vecs;
}

std::vector<std::vector<double>> brillouin_zone_path_domain<
    FERMI_SURFACE_SQUARE_2D_LATTICE>::initialize_elements() {
  std::vector<element_type> collection_k_vecs(0);

  std::vector<double> k0(2);
  std::vector<double> k1(2);

  k0[0] = M_PI;
  k0[1] = 0;
  k1[0] = 0;
  k1[1] = M_PI;

  for (int l = 0; l < INTERPOLATION_NB; l++) {
    std::vector<double> k(2, 0);

    k[0] = (1. - double(l) / double(INTERPOLATION_NB)) * k0[0] +
           double(l) / double(INTERPOLATION_NB) * k1[0];
    k[1] = (1. - double(l) / double(INTERPOLATION_NB)) * k0[1] +
           double(l) / double(INTERPOLATION_NB) * k1[1];

    collection_k_vecs.push_back(k);
  }

  collection_k_vecs.push_back(k1);

  return collection_k_vecs;
}

}  // domains
}  // phys
}  // dca
