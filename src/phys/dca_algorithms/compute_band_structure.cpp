// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements compute_band_structure.hpp.

#include "dca/phys/dca_algorithms/compute_band_structure.hpp"

namespace dca {
namespace phys {
// dca::phys::

template <>
void compute_band_structure::high_symmetry_line<1>(std::vector<std::vector<double>>& collection_k_vecs) {
  using DCA_k_cluster_type =
      domains::cluster_domain<double, 1, domains::CLUSTER, domains::MOMENTUM_SPACE,
                              domains::BRILLOUIN_ZONE>;

  int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;

  std::vector<double>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];

  for (int l = 0; l < Nb_interpolation; l++) {
    std::vector<double> k(1, 0);

    k[0] = double(l) / double(Nb_interpolation) * b0[0];

    collection_k_vecs.push_back(k);
  }
}

template <>
void compute_band_structure::high_symmetry_line<2>(std::vector<std::vector<double>>& collection_k_vecs) {
  using DCA_k_cluster_type =
      domains::cluster_domain<double, 2, domains::CLUSTER, domains::MOMENTUM_SPACE,
                              domains::BRILLOUIN_ZONE>;

  int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;

  std::vector<double>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];
  std::vector<double>& b1 = DCA_k_cluster_type::get_super_basis_vectors()[1];

  std::vector<double> k0(2);
  std::vector<double> k1(2);
  std::vector<double> k2(2);
  std::vector<double> k3(2);
  std::vector<double> k4(2);

  for (int i = 0; i < 2; i++) {
    k0[i] = 0 * b0[i] + 0. * b1[i];
    k1[i] = 1. / 2. * b0[i] + 1. / 2. * b1[i];
    k2[i] = 1. / 2. * b0[i] + 0. * b1[i];
    k3[i] = 0. * b0[i] + 1. / 2. * b1[i];
    k4[i] = 0. * b0[i] + 0. * b1[i];
  }

  for (int l = 0; l < Nb_interpolation; l++) {
    std::vector<double> k(2, 0);

    k[0] = (1. - double(l) / double(Nb_interpolation)) * k0[0] +
           double(l) / double(Nb_interpolation) * k1[0];
    k[1] = (1. - double(l) / double(Nb_interpolation)) * k0[1] +
           double(l) / double(Nb_interpolation) * k1[1];

    collection_k_vecs.push_back(k);
  }

  for (int l = 0; l < Nb_interpolation; l++) {
    std::vector<double> k(2, 0);

    k[0] = (1. - double(l) / double(Nb_interpolation)) * k1[0] +
           double(l) / double(Nb_interpolation) * k2[0];
    k[1] = (1. - double(l) / double(Nb_interpolation)) * k1[1] +
           double(l) / double(Nb_interpolation) * k2[1];

    collection_k_vecs.push_back(k);
  }

  for (int l = 0; l < Nb_interpolation; l++) {
    std::vector<double> k(2, 0);

    k[0] = (1. - double(l) / double(Nb_interpolation)) * k2[0] +
           double(l) / double(Nb_interpolation) * k3[0];
    k[1] = (1. - double(l) / double(Nb_interpolation)) * k2[1] +
           double(l) / double(Nb_interpolation) * k3[1];

    collection_k_vecs.push_back(k);
  }

  for (int l = 0; l < Nb_interpolation; l++) {
    std::vector<double> k(2, 0);

    k[0] = (1. - double(l) / double(Nb_interpolation)) * k3[0] +
           double(l) / double(Nb_interpolation) * k4[0];
    k[1] = (1. - double(l) / double(Nb_interpolation)) * k3[1] +
           double(l) / double(Nb_interpolation) * k4[1];

    collection_k_vecs.push_back(k);
  }
}

template <>
void compute_band_structure::high_symmetry_line<3>(std::vector<std::vector<double>>& collection_k_vecs) {
  using DCA_k_cluster_type =
      domains::cluster_domain<double, 3, domains::CLUSTER, domains::MOMENTUM_SPACE,
                              domains::BRILLOUIN_ZONE>;

  int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;
  std::vector<double>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];
  std::vector<double>& b1 = DCA_k_cluster_type::get_super_basis_vectors()[1];
  std::vector<double>& b2 = DCA_k_cluster_type::get_super_basis_vectors()[2];

  std::vector<double> k0(3);
  std::vector<double> k1(3);
  std::vector<double> k2(3);
  std::vector<double> k3(3);

  std::vector<double> kA(3);
  std::vector<double> kB(3);

  for (int i = 0; i < 3; i++) {
    k0[i] = 0 * b0[i] + 0. * b1[i] + 0. * b2[i];
    k1[i] = 1. / 2. * b0[i] + 1. / 2. * b1[i] + 1. / 2. * b2[i];
    k2[i] = 1. / 2. * b0[i] + 1. / 2. * b1[i] + 0. * b2[i];
    k3[i] = 1. / 2. * b0[i] + 0. * b1[i] + 0. * b2[i];

    kA[i] = -3. * b0[i] - 3. * b1[i] - 3. * b2[i];
    kB[i] = 3. * b0[i] + 3. * b1[i] + 3. * b2[i];
  }

  for (int l = 0; l < Nb_interpolation; l++) {
    std::vector<double> k(3, 0);

    k[0] = (1. - double(l) / double(Nb_interpolation)) * k0[0] +
           double(l) / double(Nb_interpolation) * k1[0];
    k[1] = (1. - double(l) / double(Nb_interpolation)) * k0[1] +
           double(l) / double(Nb_interpolation) * k1[1];
    k[2] = (1. - double(l) / double(Nb_interpolation)) * k0[2] +
           double(l) / double(Nb_interpolation) * k1[2];

    collection_k_vecs.push_back(k);
  }

  for (int l = 0; l < Nb_interpolation; l++) {
    std::vector<double> k(3, 0);

    k[0] = (1. - double(l) / double(Nb_interpolation)) * k1[0] +
           double(l) / double(Nb_interpolation) * k2[0];
    k[1] = (1. - double(l) / double(Nb_interpolation)) * k1[1] +
           double(l) / double(Nb_interpolation) * k2[1];
    k[2] = (1. - double(l) / double(Nb_interpolation)) * k1[2] +
           double(l) / double(Nb_interpolation) * k2[2];

    collection_k_vecs.push_back(k);
  }

  for (int l = 0; l < Nb_interpolation; l++) {
    std::vector<double> k(3, 0);

    k[0] = (1. - double(l) / double(Nb_interpolation)) * k2[0] +
           double(l) / double(Nb_interpolation) * k3[0];
    k[1] = (1. - double(l) / double(Nb_interpolation)) * k2[1] +
           double(l) / double(Nb_interpolation) * k3[1];
    k[2] = (1. - double(l) / double(Nb_interpolation)) * k2[2] +
           double(l) / double(Nb_interpolation) * k3[2];

    collection_k_vecs.push_back(k);
  }

  for (int l = 0; l < Nb_interpolation; l++) {
    std::vector<double> k(3, 0);

    k[0] = (1. - double(l) / double(Nb_interpolation)) * k3[0] +
           double(l) / double(Nb_interpolation) * k0[0];
    k[1] = (1. - double(l) / double(Nb_interpolation)) * k3[1] +
           double(l) / double(Nb_interpolation) * k0[1];
    k[2] = (1. - double(l) / double(Nb_interpolation)) * k3[2] +
           double(l) / double(Nb_interpolation) * k0[2];

    collection_k_vecs.push_back(k);
  }
}

template <>
void compute_band_structure::high_symmetry_line<1>(std::string& name,
                                                   std::vector<std::vector<double>>& k_vecs) {
  using DCA_k_cluster_type =
      domains::cluster_domain<double, 1, domains::CLUSTER, domains::MOMENTUM_SPACE,
                              domains::BRILLOUIN_ZONE>;

  name = "absolute";

  std::vector<double>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];

  std::vector<double> k0(1);
  std::vector<double> k1(1);

  for (int i = 0; i < 1; i++) {
    k0[i] = 0 * b0[i];        //+0.   *b1[i];
    k1[i] = 1. / 2. * b0[i];  //+1./2.*b1[i];
  }

  k_vecs.resize(0);

  k_vecs.push_back(k0);
  k_vecs.push_back(k1);
  // Should we add k0 again to close the path? (See 2D and 3D cases.)
}

template <>
void compute_band_structure::high_symmetry_line<2>(std::string& name,
                                                   std::vector<std::vector<double>>& k_vecs) {
  using DCA_k_cluster_type =
      domains::cluster_domain<double, 2, domains::CLUSTER, domains::MOMENTUM_SPACE,
                              domains::BRILLOUIN_ZONE>;

  name = "absolute";

  std::vector<double>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];
  std::vector<double>& b1 = DCA_k_cluster_type::get_super_basis_vectors()[1];

  std::vector<double> k0(2);
  std::vector<double> k1(2);
  std::vector<double> k2(2);
  std::vector<double> k3(2);

  for (int i = 0; i < 2; i++) {
    k0[i] = 0 * b0[i] + 0. * b1[i];
    k1[i] = 1. / 2. * b0[i] + 1. / 2. * b1[i];
    k2[i] = 1. / 2. * b0[i] + 0. * b1[i];
    k3[i] = 0. * b0[i] + 1. / 2. * b1[i];
  }

  k_vecs.resize(0);

  k_vecs.push_back(k0);
  k_vecs.push_back(k1);
  k_vecs.push_back(k2);
  k_vecs.push_back(k3);
  k_vecs.push_back(k0);
}

template <>
void compute_band_structure::high_symmetry_line<3>(std::string& name,
                                                   std::vector<std::vector<double>>& k_vecs) {
  using DCA_k_cluster_type =
      domains::cluster_domain<double, 3, domains::CLUSTER, domains::MOMENTUM_SPACE,
                              domains::BRILLOUIN_ZONE>;

  name = "absolute";

  std::vector<double>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];
  std::vector<double>& b1 = DCA_k_cluster_type::get_super_basis_vectors()[1];
  std::vector<double>& b2 = DCA_k_cluster_type::get_super_basis_vectors()[2];

  std::vector<double> k0(3);
  std::vector<double> k1(3);
  std::vector<double> k2(3);
  std::vector<double> k3(3);

  for (int i = 0; i < 3; i++) {
    k0[i] = 0 * b0[i] + 0. * b1[i] + 0. * b2[i];
    k1[i] = 1. / 2. * b0[i] + 1. / 2. * b1[i] + 1. / 2. * b2[i];
    k2[i] = 1. / 2. * b0[i] + 1. / 2. * b1[i] + 0. * b2[i];
    k3[i] = 1. / 2. * b0[i] + 0. * b1[i] + 0. * b2[i];
  }

  k_vecs.resize(0);

  k_vecs.push_back(k0);
  k_vecs.push_back(k1);
  k_vecs.push_back(k2);
  k_vecs.push_back(k3);
  k_vecs.push_back(k0);
}

}  // phys
}  // dca
