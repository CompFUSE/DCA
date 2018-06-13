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
// This class computes the band structure.
//
// TODO: Write out Brillouin zone cut domain, which is initialized in this class.

#ifndef DCA_PHYS_DCA_ALGORITHMS_COMPUTE_BAND_STRUCTURE_HPP
#define DCA_PHYS_DCA_ALGORITHMS_COMPUTE_BAND_STRUCTURE_HPP

#include <cassert>
#include <string>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/quantum/brillouin_zone_cut_domain.hpp"

namespace dca {
namespace phys {
// dca::phys::

class compute_band_structure {
public:
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;

  using brillouin_zone_cut_domain_type = domains::brillouin_zone_cut_domain<101>;
  using k_domain_cut_dmn_type = func::dmn_0<brillouin_zone_cut_domain_type>;
  using nu_k_cut = func::dmn_variadic<nu, k_domain_cut_dmn_type>;

  // We give it a prime number such that it is positioned on an edge or corner of a patch.
  static const int INTERPOLATION_POINTS_BAND_STRUCTURE =
      brillouin_zone_cut_domain_type::INTERPOLATION_POINTS;  // 101;

  // Computes the band structure of the non-interacting Hamiltonian H_0.
  template <typename ParametersType>
  static void execute(const ParametersType& parameters, func::function<double, nu_k_cut>& bands);

private:
  template <int lattice_dimension>
  static void construct_path(std::string coordinate_type, std::vector<std::vector<double>> path_vecs,
                             std::vector<std::vector<double>>& collection_k_vecs);

  // TODO: Pass lattice_dimension as function parameter.
  template <int lattice_dimension>
  static void high_symmetry_line(std::vector<std::vector<double>>& collection_k_vecs);

  template <int lattice_dimension>
  static void high_symmetry_line(std::string& name, std::vector<std::vector<double>>& k_vecs);
};

template <typename ParametersType>
void compute_band_structure::execute(const ParametersType& parameters,
                                     func::function<double, nu_k_cut>& band_structure) {
  std::vector<std::vector<double>> collection_k_vecs;

  std::string coordinate_type;
  std::vector<std::vector<double>> brillouin_zone_vecs;

  // Construct the path in the Brilluoin zone.
  high_symmetry_line<ParametersType::lattice_dimension>(coordinate_type, brillouin_zone_vecs);
  construct_path<ParametersType::lattice_dimension>(coordinate_type, brillouin_zone_vecs,
                                                    collection_k_vecs);

  brillouin_zone_cut_domain_type::get_size() = collection_k_vecs.size();
  brillouin_zone_cut_domain_type::get_elements() = collection_k_vecs;

  band_structure.reset();

  // Compute H(k).
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_domain_cut_dmn_type>> H_k(
      "H_k");
  ParametersType::lattice_type::initialize_H_0(parameters, H_k);

  // Compute the bands.
  dca::linalg::Vector<double, dca::linalg::CPU> L_vec(nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> H_mat(nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> V_mat(nu::dmn_size());

  for (int l = 0; l < int(collection_k_vecs.size()); l++) {
    for (int i = 0; i < nu::dmn_size(); i++)
      for (int j = 0; j < nu::dmn_size(); j++)
        H_mat(i, j) = H_k(i, j, l);

    dca::linalg::matrixop::eigensolverHermitian('N', 'U', H_mat, L_vec, V_mat);

    for (int i = 0; i < b::dmn_size(); i++)
      for (int j = 0; j < s::dmn_size(); j++)
        band_structure(i, j, l) = L_vec[2 * i + j];
  }
}

template <int lattice_dimension>
void compute_band_structure::construct_path(std::string coordinate_type,
                                            std::vector<std::vector<double>> path_vecs,
                                            std::vector<std::vector<double>>& collection_k_vecs) {
  using DCA_k_cluster_type =
      domains::cluster_domain<double, lattice_dimension, domains::CLUSTER, domains::MOMENTUM_SPACE,
                              domains::BRILLOUIN_ZONE>;

  int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;

  collection_k_vecs.resize(0);

  std::vector<double> k(lattice_dimension);

  for (int i = 0; i < int(path_vecs.size()) - 1; i++) {
    std::vector<double>& p0 = path_vecs[i];
    std::vector<double>& p1 = path_vecs[i + 1];

    assert(lattice_dimension == int(p0.size()));
    assert(lattice_dimension == int(p1.size()));

    std::vector<double> p0_tmp(lattice_dimension, 0.);
    std::vector<double> p1_tmp(lattice_dimension, 0.);

    if (coordinate_type == "relative") {
      for (int di = 0; di < lattice_dimension; di++)
        for (int dj = 0; dj < lattice_dimension; dj++)
          p0_tmp[dj] += p0[di] * DCA_k_cluster_type::get_super_basis_vectors()[di][dj];

      for (int di = 0; di < lattice_dimension; di++)
        for (int dj = 0; dj < lattice_dimension; dj++)
          p1_tmp[dj] += p1[di] * DCA_k_cluster_type::get_super_basis_vectors()[di][dj];
    }
    else {
      p0_tmp = p0;
      p1_tmp = p1;
    }

    for (int l = 0; l < Nb_interpolation; l++) {
      for (int z = 0; z < lattice_dimension; z++)
        k[z] = (1. - double(l) / double(Nb_interpolation)) * p0_tmp[z] +
               double(l) / double(Nb_interpolation) * p1_tmp[z];

      collection_k_vecs.push_back(k);
    }
  }
}

template <>
void compute_band_structure::high_symmetry_line<1>(std::vector<std::vector<double>>& collection_k_vecs);
template <>
void compute_band_structure::high_symmetry_line<2>(std::vector<std::vector<double>>& collection_k_vecs);
template <>
void compute_band_structure::high_symmetry_line<3>(std::vector<std::vector<double>>& collection_k_vecs);

template <>
void compute_band_structure::high_symmetry_line<1>(std::string& name,
                                                   std::vector<std::vector<double>>& k_vecs);
template <>
void compute_band_structure::high_symmetry_line<3>(std::string& name,
                                                   std::vector<std::vector<double>>& k_vecs);
template <>
void compute_band_structure::high_symmetry_line<3>(std::string& name,
                                                   std::vector<std::vector<double>>& k_vecs);

}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ALGORITHMS_COMPUTE_BAND_STRUCTURE_HPP
