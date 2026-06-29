// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Thomas A. Maier (maierta@ornl.gov)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Implementation of  the two-orbital model described in "Two pairing domes as Cu2+ varies to Cu3+."
// https://link.aps.org/doi/10.1103/PhysRevB.99.224515

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_THREEBAND_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_THREEBAND_LATTICE_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename SymmetryGroup>
class ThreebandHubbard {
public:
  static constexpr bool complex_g0 = false;
  static constexpr bool spin_symmetric = true;

  using LDA_point_group = domains::no_symmetry<2>;
  using DCA_point_group = SymmetryGroup;

  const static ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  const static ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS = 3;

  static const double* initializeRDCABasis();
  static const double* initializeKDCABasis();

  static const double* initializeRLDABasis();
  static const double* initializeKLDABasis();

  static std::vector<int> flavors();
  static std::vector<std::vector<double>> aVectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> orbitalPermutations();

  static int transformationSignOfR(int b1, int b2, int s);
  static int transformationSignOfK(int b1, int b2, int s);

  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initializeHInteraction(
      func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
      const parameters_type& parameters);

  template <class domain>
  static void initializeHSymmetry(func::function<int, domain>& H_symmetry);

  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initializeH0(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
      func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);

  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initializeH0WithQ(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
      func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0, typename KDmn::element_type q);
};

template <typename PointGroupType>
int ThreebandHubbard<PointGroupType>::transformationSignOfR(int b1, int b2, int s) {
  if (!std::is_same<PointGroupType, domains::D4>::value)
    return 1;

  if (b1 == b2)
    return 1;
  else if (b1 != 0 && b2 != 0) {
    if (s == 0 || s == 6)
      return 1;
    else
      return 0;
  }
  else if ((b1 != b2 && b1 == 0) || (b1 != b2 && b2 == 0)) {
    if (s == 0)
      return 1;
    else if (s == 6)
      return -1;
    else
      return 0;
  }

  return s == 0;
}

template <typename PointGroupType>
int ThreebandHubbard<PointGroupType>::transformationSignOfK(int b1, int b2, int s) {
  if (!std::is_same<PointGroupType, domains::D4>::value)
    return 1;

  if ((b1 == b2) || (b1 != 0 && b2 != 0))
    return 1;
  else
    return (s == 0 || s == 2 || s == 5 || s == 7) ? 1 : -1;
}

template <typename PointGroupType>
const double* ThreebandHubbard<PointGroupType>::initializeRDCABasis() {
  static std::array<double, 4> basis{1, 0, 0, 1};
  return basis.data();
}

template <typename PointGroupType>
const double* ThreebandHubbard<PointGroupType>::initializeKDCABasis() {
  static std::array<double, 4> basis{2 * M_PI, 0, 0, 2 * M_PI};
  return basis.data();
}

template <typename PointGroupType>
const double* ThreebandHubbard<PointGroupType>::initializeRLDABasis() {
  static std::array<double, 4> basis{1, 0, 0, 1};
  return basis.data();
}

template <typename PointGroupType>
const double* ThreebandHubbard<PointGroupType>::initializeKLDABasis() {
  static std::array<double, 4> basis{2 * M_PI, 0, 0, 2 * M_PI};
  return basis.data();
}

template <typename PointGroupType>
std::vector<int> ThreebandHubbard<PointGroupType>::flavors() {
  return {0, 1, 1};
}

template <typename PointGroupType>
std::vector<std::vector<double>> ThreebandHubbard<PointGroupType>::aVectors() {
  return {{0, 0}, {0.5, 0}, {0, 0.5}};
}

template <typename PointGroupType>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>
ThreebandHubbard<PointGroupType>::orbitalPermutations() {
  return {};
}

template <typename PointGroupType>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void ThreebandHubbard<PointGroupType>::initializeHInteraction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Three-band Hubbard model has three bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double U_dd = parameters.get_U_dd();
  const double U_pp = parameters.get_U_pp();
  const double V_pd = parameters.get_V_pd();
  const double V_pp = parameters.get_V_pp();

  H_interaction = 0.;

  const auto& basis = RDmn::parameter_type::get_basis_vectors();
  const auto& super_basis = RDmn::parameter_type::get_super_basis_vectors();
  const auto& elements = RDmn::parameter_type::get_elements();

  if (basis.size() != DIMENSION)
    throw std::logic_error("Unexpected lattice basis dimension for three-band Hubbard model.");

  auto index_of = [&](std::vector<double> displacement) {
    displacement = domains::cluster_operations::translate_inside_cluster(displacement, super_basis);
    return domains::cluster_operations::index(displacement, elements, domains::BRILLOUIN_ZONE);
  };

  auto push_unique = [](std::vector<int>& indices, int index) {
    if (std::find(indices.begin(), indices.end(), index) == indices.end())
      indices.push_back(index);
  };

  const int plus_x = index_of(basis[0]);
  const int plus_y = index_of(basis[1]);
  const int minus_x = RDmn::parameter_type::subtract(plus_x, origin);
  const int minus_y = RDmn::parameter_type::subtract(plus_y, origin);

  std::vector<double> plus_x_minus_y(DIMENSION, 0.);
  std::vector<double> minus_x_plus_y(DIMENSION, 0.);
  for (int d = 0; d < DIMENSION; ++d) {
    plus_x_minus_y[d] = basis[0][d] - basis[1][d];
    minus_x_plus_y[d] = basis[1][d] - basis[0][d];
  }

  const int delta_plus_x_minus_y = index_of(plus_x_minus_y);
  const int delta_minus_x_plus_y = index_of(minus_x_plus_y);

  std::vector<int> d_to_px{origin};
  std::vector<int> px_to_d{origin};
  std::vector<int> d_to_py{origin};
  std::vector<int> py_to_d{origin};
  std::vector<int> px_to_py{origin};
  std::vector<int> py_to_px{origin};

  push_unique(d_to_px, plus_x);
  push_unique(px_to_d, minus_x);
  push_unique(d_to_py, plus_y);
  push_unique(py_to_d, minus_y);

  push_unique(px_to_py, minus_x);
  push_unique(px_to_py, plus_y);
  push_unique(px_to_py, delta_minus_x_plus_y);

  push_unique(py_to_px, plus_x);
  push_unique(py_to_px, minus_y);
  push_unique(py_to_px, delta_plus_x_minus_y);

  for (int i = 0; i < BANDS; i++) {
    for (int s1 = 0; s1 < 2; s1++) {
      for (int j = 0; j < BANDS; j++) {
        for (int s2 = 0; s2 < 2; s2++) {
          if (i == 0 && j == 0 && s1 != s2)
            H_interaction(i, s1, j, s2, origin) = U_dd;

          if (i == j && i != 0 && s1 != s2)
            H_interaction(i, s1, j, s2, origin) = U_pp;
        }
      }
    }
  }

  for (int s1 = 0; s1 < 2; s1++) {
    for (int s2 = 0; s2 < 2; s2++) {
      for (const int delta_r : d_to_px)
        H_interaction(0, s1, 1, s2, delta_r) = V_pd;
      for (const int delta_r : d_to_py)
        H_interaction(0, s1, 2, s2, delta_r) = V_pd;
      for (const int delta_r : px_to_d)
        H_interaction(1, s1, 0, s2, delta_r) = V_pd;
      for (const int delta_r : py_to_d)
        H_interaction(2, s1, 0, s2, delta_r) = V_pd;
    }
  }

  for (int s1 = 0; s1 < 2; s1++) {
    for (int s2 = 0; s2 < 2; s2++) {
      for (const int delta_r : px_to_py)
        H_interaction(1, s1, 2, s2, delta_r) = V_pp;
      for (const int delta_r : py_to_px)
        H_interaction(2, s1, 1, s2, delta_r) = V_pp;
    }
  }
}

template <typename PointGroupType>
template <class domain>
void ThreebandHubbard<PointGroupType>::initializeHSymmetry(
    func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;
}

template <typename PointGroupType>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void ThreebandHubbard<PointGroupType>::initializeH0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  typename KDmn::element_type default_q(ParametersType::lattice_dimension);
  initializeH0WithQ(parameters, H_0, default_q);
}

template <typename PointGroupType>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void ThreebandHubbard<PointGroupType>::initializeH0WithQ(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0, typename KDmn::element_type q) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Three-band Hubbard model has three bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  std::vector<typename KDmn::element_type> k_vecs(KDmn::get_elements());

  for (auto& k_elem : k_vecs)
    std::transform(k_elem.cbegin(), k_elem.cend(), q.begin(), k_elem.begin(),
                   [](auto& k, auto q_elem) { return k + q_elem; });

  const auto t_pd = parameters.get_t_pd();
  const auto t_pp = parameters.get_t_pp();
  const auto ep_d = parameters.get_ep_d();
  const auto ep_p = parameters.get_ep_p();
  const auto U_dd = parameters.get_U_dd();
  const auto U_pp = parameters.get_U_pp();
  const auto V_pd = parameters.get_V_pd();
  const auto V_pp = parameters.get_V_pp();

  H_0 = ScalarType(0);

  const ScalarType I(0, 1.);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];

    const auto valdpx = -2. * I * t_pd * std::sin(k[0] / 2.);
    const auto valdpy = 2. * I * t_pd * std::sin(k[1] / 2.);
    const auto valpxpy = 4. * t_pp * std::sin(k[0] / 2.) * std::sin(k[1] / 2.);

    for (int s = 0; s < 2; s++) {
      H_0(0, s, 0, s, k_ind) = ep_d + U_dd / 2.0 + 4.0 * V_pd;
      H_0(1, s, 1, s, k_ind) = ep_p + U_pp / 2.0 + 4.0 * V_pp + 2.0 * V_pd;
      H_0(2, s, 2, s, k_ind) = ep_p + U_pp / 2.0 + 4.0 * V_pp + 2.0 * V_pd;

      H_0(0, s, 1, s, k_ind) = valdpx;
      H_0(1, s, 0, s, k_ind) = -valdpx;

      H_0(0, s, 2, s, k_ind) = valdpy;
      H_0(2, s, 0, s, k_ind) = -valdpy;

      H_0(1, s, 2, s, k_ind) = valpxpy;
      H_0(2, s, 1, s, k_ind) = valpxpy;
    }
  }
}

}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_THREEBAND_LATTICE_HPP
