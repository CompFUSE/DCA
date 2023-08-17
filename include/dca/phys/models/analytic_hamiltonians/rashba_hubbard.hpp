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
// Square lattice.

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_RASHBA_HUBBARD_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_RASHBA_HUBBARD_HPP

#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/util.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename PointGroup>
class RashbaHubbard {
public:
  static constexpr bool complex_g0 = true;
  static constexpr bool spin_symmetric = false;

  using LDA_point_group = domains::no_symmetry<2>;
  using DCA_point_group = PointGroup;

  const static int DIMENSION = 2;

  const static int SPINS = 2;
  // The model is singleband, but up and down spins are stored in the same sector.
  const static int BANDS = 2;

  static const double* initializeRDCABasis();
  static const double* initializeRLDABasis();

  constexpr static int transformationSignOfR(int, int, int) {
    return 1;
  }
  constexpr static int transformationSignOfK(int, int, int) {
    return 1;
  }
  // constexpr static int H0UpDownFinite(int k) {
  //   // return 1;
  //   return H0UpDownFiniteK(k);
  // }

  static std::vector<int> flavors();
  static std::vector<int> spins();

  static std::vector<std::vector<double>> aVectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> orbitalPermutations();

  // Initializes the interaction part of the real space Hubbard Hamiltonian.
  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initializeHInteraction(
      func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
      const parameters_type& parameters);

  template <class domain>
  static void initializeHSymmetry(func::function<int, domain>& H_symmetry);

  // Initializes the tight-binding (non-interacting) part of the momentum space Hamiltonian.
  // Preconditions: The elements of KDmn are two-dimensional (access through index 0 and 1).
  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initializeH0(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);
  template <typename FNEW, typename KVECS>
  static void timeOrFrequencySymmetrySpecial(dca::phys::domains::CLUSTER_REPRESENTATION cr,
                                             FNEW& function, KVECS& k_vecs, int c_ind, int w_ind,
                                             int w_0);

  static void clusterSymmetrySpecial(int b0, int b1, int k_ind, int& k_new, int& b0_new,
                                     int& b1_new, double& sign);
};

template <typename PointGroup>
const double* RashbaHubbard<PointGroup>::initializeRDCABasis() {
  static const std::array<double, 4> r_base{1, 0, 0, 1};
  return r_base.data();
}

template <typename PointGroup>
const double* RashbaHubbard<PointGroup>::initializeRLDABasis() {
  static const std::array<double, 4> r_base{1, 0, 0, 1};
  return r_base.data();
}

template <typename PointGroup>
std::vector<int> RashbaHubbard<PointGroup>::spins() {
  static std::vector<int> spins(spins);

  for (int i = 0; i < SPINS; i++)
    spins[i] = i;

  return spins;
}
  
template <typename PointGroup>
std::vector<int> RashbaHubbard<PointGroup>::flavors() {
  static std::vector<int> flavors(BANDS);

  for (int i = 0; i < BANDS; i++)
    flavors[i] = i;

  return flavors;
}

template <typename PointGroup>
std::vector<std::vector<double>> RashbaHubbard<PointGroup>::aVectors() {
  return {{0, 0}, {0, 0}};
}

template <typename PointGroup>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> RashbaHubbard<
    PointGroup>::orbitalPermutations() {
  return {};
}

template <typename PointGroup>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void RashbaHubbard<PointGroup>::initializeHInteraction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Rashba lattice has 2 bands.");
  // if (SpinDmn::dmn_size() != 2)
  //   throw std::logic_error("Spin domain size must be 2.");

  // Get the index of the origin (0,0).
  const int origin = RDmn::parameter_type::origin_index();

  // Set all elements to zero.
  H_interaction = 0.;

  // On-site interaction, store up-down interaction in the first sector.
  const double U = parameters.get_U();
  H_interaction(0, 0, 1, 0, origin) = U;
  H_interaction(1, 0, 0, 0, origin) = U;
}

template <typename PointGroup>
template <class domain>
void RashbaHubbard<PointGroup>::initializeHSymmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  // H_symmetries(0, 0, 0, 0) = 0; // at b0, G of spin 0 or 1 has the same values.
  // H_symmetries(1, 0, 1, 0) = 0;
  // H_symmetries(0, 1, 0, 1) = 0; // at b0, G of spin 0 or 1 has the same values.
  // H_symmetries(1, 0, 1, 0) = 0;
}

template <typename PointGroup>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void RashbaHubbard<PointGroup>::initializeH0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Square lattice has one band.");
  // if (SpinDmn::dmn_size() != 2)
  //   throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto t = parameters.get_t();
  const auto h = parameters.get_h();
  const auto lambda = parameters.get_lambda();

  H_0 = 0;
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> m(2);
  constexpr ScalarType i{0, 1};

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];

    // kinetic term
    m(0, 0) = m(1, 1) = -2. * t * (std::cos(k[0]) + std::cos(k[1]));

    // Note: spin space is {e_DN, e_UP}
    // Zeeman field
    m(0, 0) += h;
    m(1, 1) += -h;

    // Pauli matrix sigma_x
    m(0, 1) = m(1, 0) = 2 * lambda * (-std::sin(k[1]));

    // Pauli matrix sigma_y
    const auto val = 2 * lambda * std::sin(k[0]);
    m(0, 1) += +i * val;
    m(1, 0) += -i * val;

    for (int s1 = 0; s1 < 2; ++s1)
      for (int s2 = 0; s2 < 2; ++s2)
        H_0(s1, 0, s2, 0, k_ind) = m(s1, s2);
  }
}

// This breaks single band models symmetrized over spin and probably produces something completely
// wrong in other cases For Rashba model: Set inter-orbital (spin-up/down) component to zero when
// sin(kx)=0 & sin(ky)=0, i.e. when inter-orbital (inter-spin) Hamiltonian is zero
template <typename PointGroup>
template <typename FNEW, typename KVECS>
void RashbaHubbard<PointGroup>::timeOrFrequencySymmetrySpecial(
    dca::phys::domains::CLUSTER_REPRESENTATION cr, FNEW& f_new, KVECS& k_vecs, int c_ind, int w_ind,
    int w_0) {
  if (cr == domains::MOMENTUM_SPACE) {
    const auto& k = k_vecs[c_ind];

    if (abs(std::sin(k[0])) < 1.0e-4 && abs(std::sin(k[1])) < 1.0e-4) {
      // std::cout << "Setting off-diag comp. to zero\n";

      f_new(0, 1, c_ind, w_ind) = 0.0;
      f_new(1, 0, c_ind, w_ind) = 0.0;
      f_new(0, 1, c_ind, w_0 - w_ind) = 0.0;
      f_new(1, 0, c_ind, w_0 - w_ind) = 0.0;
    }
  }
}

template <typename PointGroup>
void RashbaHubbard<PointGroup>::clusterSymmetrySpecial(int b0, int b1, int k_ind, int& k_new,
                                                       int& b0_new, int& b1_new, double& sign) {
  if (b0 != b1) {  // For Rashba model, the up-down elements transform like ix + y
    // const auto& k_vecs = func::dmn_0<domains::cluster_domain<scalar_type, D, N,
    // domains::MOMENTUM_SPACE, S>>::get_elements(); const auto& k1 = k_vecs[k_ind]; const
    // auto& k2 = k_vecs[k_new];
    k_new = k_ind;
    b0_new = b0;
    b1_new = b1;
    sign = 1;
  }
}

// template <typename PointGroup>
// template <typename KDmn>
// void RashbaHubbard<PointGroup>::H0UpDownFiniteK(int k) {

// }

}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_RASHBA_HUBBARD_HPP
