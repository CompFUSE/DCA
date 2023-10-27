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

#include "dca/config/haves_defines.hpp"
#include "dca/platform/gpu_definitions.h"
#include "dca/platform/dca_gpu.h"

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

/** class to calculate band structure
 *  requirements have changed from when it was written.
 *  Now it could be just templated functions with access to a type alias utility struct
 */
template <class PARAMETERS>
class compute_band_structure {
public:
  using Parameters = PARAMETERS;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;

  using Real = typename Parameters::Real;

  using brillouin_zone_cut_domain_type = domains::brillouin_zone_cut_domain<101>;
  using k_domain_cut_dmn_type = func::dmn_0<brillouin_zone_cut_domain_type>;
  using nu_k_cut = func::dmn_variadic<nu, k_domain_cut_dmn_type>;

  // We give it a prime number such that it is positioned on an edge or corner of a patch.
  static const int INTERPOLATION_POINTS_BAND_STRUCTURE =
      brillouin_zone_cut_domain_type::INTERPOLATION_POINTS;  // 101;

  // Computes the band structure of the non-interacting Hamiltonian H_0.
  static void execute(const Parameters& parameters, func::function<Real, nu_k_cut>& bands);

private:
  template <int lattice_dimension>
  static void construct_path(std::string coordinate_type, std::vector<std::vector<Real>> path_vecs,
                             std::vector<std::vector<Real>>& collection_k_vecs);

  // TODO: Pass lattice_dimension as function parameter.
  template <int lattice_dimension>
  static void high_symmetry_line(std::vector<std::vector<Real>>& collection_k_vecs);

  template <int lattice_dimension>
  static void high_symmetry_line(std::string& name, std::vector<std::vector<Real>>& k_vecs);
};

template <typename PARAMETERS>
void compute_band_structure<PARAMETERS>::execute(const Parameters& parameters,
                                                 func::function<Real, nu_k_cut>& band_structure) {
  if constexpr (std::is_same_v<double, Real>) {
    std::vector<std::vector<double>> collection_k_vecs;

    std::string coordinate_type;
    std::vector<std::vector<double>> brillouin_zone_vecs;

    // Construct the path in the Brilluoin zone.
    high_symmetry_line<Parameters::lattice_dimension>(coordinate_type, brillouin_zone_vecs);
    construct_path<Parameters::lattice_dimension>(coordinate_type, brillouin_zone_vecs,
                                                      collection_k_vecs);

    brillouin_zone_cut_domain_type::get_size() = collection_k_vecs.size();
    brillouin_zone_cut_domain_type::get_elements() = collection_k_vecs;

    band_structure.reset();

    // Compute H(k).
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_domain_cut_dmn_type>> H_k(
        "H_k");
    Parameters::lattice_type::initializeH0(parameters, H_k);

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
  else {
    std::vector<std::vector<Real>> collection_k_vecs;

    std::string coordinate_type;
    std::vector<std::vector<Real>> brillouin_zone_vecs;

    // Construct the path in the Brilluoin zone.
    high_symmetry_line<Parameters::lattice_dimension>(coordinate_type, brillouin_zone_vecs);
    construct_path<Parameters::lattice_dimension>(coordinate_type, brillouin_zone_vecs,
                                                  collection_k_vecs);

    // brillouin_zone_cut_domain_type is global state. yuck.
    brillouin_zone_cut_domain_type::get_size() = collection_k_vecs.size();
    // so we are going to need to widen the kvectors I guess.
    auto k2Double = [](auto& kvec) -> std::vector<double> {
      std::vector<double> k_converted(kvec.size());
      std::transform(kvec.begin(), kvec.end(), k_converted.begin(),
                     [](auto& val) -> typename decltype(k_converted)::value_type {
                       return static_cast<typename decltype(k_converted)::value_type>(val);
                     });
      return k_converted;
    };

    std::vector<std::vector<double>> coll_k_vecs_dbl(collection_k_vecs.size());

    for (int i = 0; i < collection_k_vecs.size(); ++i) {
      coll_k_vecs_dbl[i] = k2Double(collection_k_vecs[i]);
    }

    brillouin_zone_cut_domain_type::get_elements() = coll_k_vecs_dbl;

    band_structure.reset();

    // Compute H(k).
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_domain_cut_dmn_type>> H_k(
        "H_k");
    Parameters::lattice_type::initializeH0(parameters, H_k);

    // Compute the bands.
    dca::linalg::Vector<Real, dca::linalg::CPU> L_vec(nu::dmn_size());
    dca::linalg::Matrix<std::complex<Real>, dca::linalg::CPU> H_mat(nu::dmn_size());
    dca::linalg::Matrix<std::complex<Real>, dca::linalg::CPU> V_mat(nu::dmn_size());

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
}

template <class PARAMETERS>
template <int lattice_dimension>
void compute_band_structure<PARAMETERS>::construct_path(
    std::string coordinate_type, std::vector<std::vector<Real>> path_vecs,
    std::vector<std::vector<Real>>& collection_k_vecs) {
  using DCA_k_cluster_type =
      domains::cluster_domain<double, lattice_dimension, domains::CLUSTER, domains::MOMENTUM_SPACE,
                              domains::BRILLOUIN_ZONE>;

  int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;

  collection_k_vecs.resize(0);

  std::vector<Real> k(lattice_dimension);

  for (int i = 0; i < int(path_vecs.size()) - 1; i++) {
    std::vector<Real>& p0 = path_vecs[i];
    std::vector<Real>& p1 = path_vecs[i + 1];

    assert(lattice_dimension == int(p0.size()));
    assert(lattice_dimension == int(p1.size()));

    std::vector<Real> p0_tmp(lattice_dimension, 0.);
    std::vector<Real> p1_tmp(lattice_dimension, 0.);

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
        k[z] = (1. - static_cast<Real>(l) / static_cast<Real>(Nb_interpolation)) * p0_tmp[z] +
               static_cast<Real>(l) / static_cast<Real>(Nb_interpolation) * p1_tmp[z];

      collection_k_vecs.push_back(k);
    }
  }
}

// template <>
// void compute_band_structure::high_symmetry_line<1>(std::vector<std::vector<double>>& collection_k_vecs);
// template <>
// void compute_band_structure::high_symmetry_line<2>(std::vector<std::vector<double>>& collection_k_vecs);
// template <>
// void compute_band_structure::high_symmetry_line<3>(std::vector<std::vector<double>>& collection_k_vecs);

// template <>
// void compute_band_structure::high_symmetry_line<1>(std::string& name,
//                                                    std::vector<std::vector<double>>& k_vecs);
// template <>
// void compute_band_structure::high_symmetry_line<3>(std::string& name,
//                                                    std::vector<std::vector<double>>& k_vecs);
// template <>
// void compute_band_structure::high_symmetry_line<3>(std::string& name,
//                                                    std::vector<std::vector<double>>& k_vecs);

template <class PARAMETERS>
template <int lattice_dimension>
void compute_band_structure<PARAMETERS>::high_symmetry_line(
    std::vector<std::vector<Real>>& collection_k_vecs) {
  if constexpr (lattice_dimension == 1) {
    using DCA_k_cluster_type =
        domains::cluster_domain<Real, 1, domains::CLUSTER, domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;

    int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;

    std::vector<Real>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];

    for (int l = 0; l < Nb_interpolation; l++) {
      std::vector<Real> k(1, 0);

      k[0] = static_cast<Real>(l) / static_cast<Real>(Nb_interpolation) * b0[0];

      collection_k_vecs.push_back(k);
    }
  }
  else if constexpr (lattice_dimension == 2) {
    using DCA_k_cluster_type =
        domains::cluster_domain<Real, 2, domains::CLUSTER, domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;

    int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;

    std::vector<Real>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];
    std::vector<Real>& b1 = DCA_k_cluster_type::get_super_basis_vectors()[1];

    std::vector<Real> k0(2);
    std::vector<Real> k1(2);
    std::vector<Real> k2(2);
    std::vector<Real> k3(2);
    std::vector<Real> k4(2);

    for (int i = 0; i < 2; i++) {
      k0[i] = 0 * b0[i] + 0. * b1[i];
      k1[i] = 1. / 2. * b0[i] + 1. / 2. * b1[i];
      k2[i] = 1. / 2. * b0[i] + 0. * b1[i];
      k3[i] = 0. * b0[i] + 1. / 2. * b1[i];
      k4[i] = 0. * b0[i] + 0. * b1[i];
    }

    for (int l = 0; l < Nb_interpolation; l++) {
      std::vector<Real> k(2, 0);

      k[0] =
          (1. - Real(l) / Real(Nb_interpolation)) * k0[0] + Real(l) / Real(Nb_interpolation) * k1[0];
      k[1] =
          (1. - Real(l) / Real(Nb_interpolation)) * k0[1] + Real(l) / Real(Nb_interpolation) * k1[1];

      collection_k_vecs.push_back(k);
    }

    for (int l = 0; l < Nb_interpolation; l++) {
      std::vector<Real> k(2, 0);

      k[0] =
          (1. - Real(l) / Real(Nb_interpolation)) * k1[0] + Real(l) / Real(Nb_interpolation) * k2[0];
      k[1] =
          (1. - Real(l) / Real(Nb_interpolation)) * k1[1] + Real(l) / Real(Nb_interpolation) * k2[1];

      collection_k_vecs.push_back(k);
    }

    for (int l = 0; l < Nb_interpolation; l++) {
      std::vector<Real> k(2, 0);

      k[0] =
          (1. - Real(l) / Real(Nb_interpolation)) * k2[0] + Real(l) / Real(Nb_interpolation) * k3[0];
      k[1] =
          (1. - Real(l) / Real(Nb_interpolation)) * k2[1] + Real(l) / Real(Nb_interpolation) * k3[1];

      collection_k_vecs.push_back(k);
    }

    for (int l = 0; l < Nb_interpolation; l++) {
      std::vector<Real> k(2, 0);

      k[0] =
          (1. - Real(l) / Real(Nb_interpolation)) * k3[0] + Real(l) / Real(Nb_interpolation) * k4[0];
      k[1] =
          (1. - Real(l) / Real(Nb_interpolation)) * k3[1] + Real(l) / Real(Nb_interpolation) * k4[1];

      collection_k_vecs.push_back(k);
    }
  }
  else if constexpr (lattice_dimension == 3) {
    using DCA_k_cluster_type =
        domains::cluster_domain<Real, 3, domains::CLUSTER, domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;

    int Nb_interpolation = INTERPOLATION_POINTS_BAND_STRUCTURE;
    std::vector<Real>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];
    std::vector<Real>& b1 = DCA_k_cluster_type::get_super_basis_vectors()[1];
    std::vector<Real>& b2 = DCA_k_cluster_type::get_super_basis_vectors()[2];

    std::vector<Real> k0(3);
    std::vector<Real> k1(3);
    std::vector<Real> k2(3);
    std::vector<Real> k3(3);

    std::vector<Real> kA(3);
    std::vector<Real> kB(3);

    for (int i = 0; i < 3; i++) {
      k0[i] = 0 * b0[i] + 0. * b1[i] + 0. * b2[i];
      k1[i] = 1. / 2. * b0[i] + 1. / 2. * b1[i] + 1. / 2. * b2[i];
      k2[i] = 1. / 2. * b0[i] + 1. / 2. * b1[i] + 0. * b2[i];
      k3[i] = 1. / 2. * b0[i] + 0. * b1[i] + 0. * b2[i];

      kA[i] = -3. * b0[i] - 3. * b1[i] - 3. * b2[i];
      kB[i] = 3. * b0[i] + 3. * b1[i] + 3. * b2[i];
    }

    for (int l = 0; l < Nb_interpolation; l++) {
      std::vector<Real> k(3, 0);

      k[0] =
          (1. - Real(l) / Real(Nb_interpolation)) * k0[0] + Real(l) / Real(Nb_interpolation) * k1[0];
      k[1] =
          (1. - Real(l) / Real(Nb_interpolation)) * k0[1] + Real(l) / Real(Nb_interpolation) * k1[1];
      k[2] =
          (1. - Real(l) / Real(Nb_interpolation)) * k0[2] + Real(l) / Real(Nb_interpolation) * k1[2];

      collection_k_vecs.push_back(k);
    }

    for (int l = 0; l < Nb_interpolation; l++) {
      std::vector<Real> k(3, 0);

      k[0] =
          (1. - Real(l) / Real(Nb_interpolation)) * k1[0] + Real(l) / Real(Nb_interpolation) * k2[0];
      k[1] =
          (1. - Real(l) / Real(Nb_interpolation)) * k1[1] + Real(l) / Real(Nb_interpolation) * k2[1];
      k[2] =
          (1. - Real(l) / Real(Nb_interpolation)) * k1[2] + Real(l) / Real(Nb_interpolation) * k2[2];

      collection_k_vecs.push_back(k);
    }

    for (int l = 0; l < Nb_interpolation; l++) {
      std::vector<Real> k(3, 0);

      k[0] =
          (1. - Real(l) / Real(Nb_interpolation)) * k2[0] + Real(l) / Real(Nb_interpolation) * k3[0];
      k[1] =
          (1. - Real(l) / Real(Nb_interpolation)) * k2[1] + Real(l) / Real(Nb_interpolation) * k3[1];
      k[2] =
          (1. - Real(l) / Real(Nb_interpolation)) * k2[2] + Real(l) / Real(Nb_interpolation) * k3[2];

      collection_k_vecs.push_back(k);
    }

    for (int l = 0; l < Nb_interpolation; l++) {
      std::vector<Real> k(3, 0);

      k[0] =
          (1. - Real(l) / Real(Nb_interpolation)) * k3[0] + Real(l) / Real(Nb_interpolation) * k0[0];
      k[1] =
          (1. - Real(l) / Real(Nb_interpolation)) * k3[1] + Real(l) / Real(Nb_interpolation) * k0[1];
      k[2] =
          (1. - Real(l) / Real(Nb_interpolation)) * k3[2] + Real(l) / Real(Nb_interpolation) * k0[2];

      collection_k_vecs.push_back(k);
    }
  }
}

template <class PARAMETERS>
template <int lattice_dimension>
void compute_band_structure<PARAMETERS>::high_symmetry_line(std::string& name,
                                                            std::vector<std::vector<Real>>& k_vecs) {
  if constexpr (lattice_dimension == 1) {
    using DCA_k_cluster_type =
        domains::cluster_domain<Real, 1, domains::CLUSTER, domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;

    name = "absolute";

    std::vector<Real>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];

    std::vector<Real> k0(1);
    std::vector<Real> k1(1);

    for (int i = 0; i < 1; i++) {
      k0[i] = 0 * b0[i];        //+0.   *b1[i];
      k1[i] = 1. / 2. * b0[i];  //+1./2.*b1[i];
    }

    k_vecs.resize(0);

    k_vecs.push_back(k0);
    k_vecs.push_back(k1);
    // Should we add k0 again to close the path? (See 2D and 3D cases.)
  }
  else if constexpr (lattice_dimension == 2) {
    using DCA_k_cluster_type =
        domains::cluster_domain<Real, 2, domains::CLUSTER, domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;

    name = "absolute";

    std::vector<Real>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];
    std::vector<Real>& b1 = DCA_k_cluster_type::get_super_basis_vectors()[1];

    std::vector<Real> k0(2);
    std::vector<Real> k1(2);
    std::vector<Real> k2(2);
    std::vector<Real> k3(2);

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
  else if constexpr (lattice_dimension == 3) {
    using DCA_k_cluster_type =
        domains::cluster_domain<Real, 3, domains::CLUSTER, domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;

    name = "absolute";

    std::vector<Real>& b0 = DCA_k_cluster_type::get_super_basis_vectors()[0];
    std::vector<Real>& b1 = DCA_k_cluster_type::get_super_basis_vectors()[1];
    std::vector<Real>& b2 = DCA_k_cluster_type::get_super_basis_vectors()[2];

    std::vector<Real> k0(3);
    std::vector<Real> k1(3);
    std::vector<Real> k2(3);
    std::vector<Real> k3(3);

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
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_ALGORITHMS_COMPUTE_BAND_STRUCTURE_HPP
