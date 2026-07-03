// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Tyler Sax (tylersax@gmail.com)
//
// Derives a 2D model's point group -- the built-in 2D holohedry pool filtered by the lattice geometry
// and by H0 invariance -- and compares it against the group realized from the model's declared
// DCA_point_group. Production symmetrization continues to use the DECLARED group:
// deriveAndComparePointGroup restores the declared symmetry state before returning and only
// reports divergences. The reports build a per-model map of models whose declaration disagrees with
// the symmetry encoded in their own Hamiltonian, which drives later declaration-free upgrades.

#ifndef DCA_PHYS_DCA_STEP_SYMMETRIZATION_DERIVE_POINT_GROUP_HPP
#define DCA_PHYS_DCA_STEP_SYMMETRIZATION_DERIVE_POINT_GROUP_HPP

#include <cmath>
#include <complex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/dca_step/symmetrization/solve_orbital_op_signs.hpp"
#include "dca/phys/domains/cluster/cluster_symmetry.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/holohedries_2d.hpp"
#include "dca/phys/domains/cluster/symmetrization_algorithms/search_maximal_symmetry_group.hpp"
#include "dca/phys/domains/cluster/symmetrization_algorithms/set_symmetry_matrices.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/quantum/point_group_symmetry_element.hpp"

namespace dca {
namespace phys {
// dca::phys::

namespace detail {
// dca::phys::detail::

// Geometric equality of two symmetry elements (same linear part and fractional translation), at the
// tolerance the symmetry search uses to deduplicate operations.
inline bool isSameSymmetryElement(const domains::point_group_symmetry_element& a,
                                  const domains::point_group_symmetry_element& b) {
  constexpr double tol = 1e-6;
  if (a.DIMENSION != b.DIMENSION)
    return false;
  double diff = 0.;
  for (int i = 0; i < a.DIMENSION * a.DIMENSION; ++i)
    diff = std::max(diff, std::abs(a.O[i] - b.O[i]));
  for (int i = 0; i < a.DIMENSION; ++i)
    diff = std::max(diff, std::abs(a.t[i] - b.t[i]));
  return diff < tol;
}

inline bool containsSymmetryElement(const std::vector<domains::point_group_symmetry_element>& set,
                                    const domains::point_group_symmetry_element& op) {
  for (const auto& element : set)
    if (isSameSymmetryElement(element, op))
      return true;
  return false;
}

// Detects the initializeH0 lattice interface: the legacy lattices cannot feed the
// H0-invariance gate, so derivation skips them.
template <typename Lattice, typename Parameters, typename H0Type, typename = void>
struct LatticeHasInitializeH0 : std::false_type {};
template <typename Lattice, typename Parameters, typename H0Type>
struct LatticeHasInitializeH0<Lattice, Parameters, H0Type,
                              std::void_t<decltype(Lattice::initializeH0(
                                  std::declval<const Parameters&>(), std::declval<H0Type&>()))>>
    : std::true_type {};

// Makes PointGroup the live group of one cluster family, overwriting the previously realized group:
// clear the operation lists (the search only appends), re-search under the same geometric filter,
// resize the per-operation functions to the new count, and rebuild the symmetry tables. Assumes the
// family was already populated once by cluster_domain_symmetry_initializer.
template <typename KCluster, typename PointGroup>
void installPointGroup() {
  using CS_k = domains::cluster_symmetry<KCluster>;
  using CS_r = domains::cluster_symmetry<typename CS_k::dual_type>;
  using Family = typename CS_k::cluster_family_type;

  CS_k::sym_unit_cell_dmn_t::parameter_type::get_elements().clear();
  CS_k::sym_super_cell_dmn_t::parameter_type::get_elements().clear();
  domains::search_maximal_symmetry_group<Family, PointGroup, domains::UNIT_CELL>::execute();
  domains::search_maximal_symmetry_group<Family, PointGroup, domains::SUPER_CELL>::execute();
  CS_k::get_symmetry_matrix().reset();
  CS_r::get_symmetry_matrix().reset();
  CS_k::get_mapped_point().reset();
  CS_k::get_fold_phase().reset();
  CS_k::get_orbital_op().reset();
  domains::set_symmetry_matrices<Family>::execute();
}

}  // namespace detail
// dca::phys::

// Derives the point group of a 2D model for one cluster family and returns a human-readable report
// of how it diverges from the declared group -- empty when they agree. Must be called right after
// cluster_domain_symmetry_initializer<RClusterDmn, DCA_point_group>::execute() for the same family;
// the declared symmetry state is restored before returning, so production behavior is unchanged.
// 3D models (no 3D pool yet) and legacy lattices (no initializeH0) skip derivation entirely.
template <typename RClusterDmn, typename Model, typename Parameters>
std::string deriveAndComparePointGroup(const Parameters& parameters) {
  using Lattice = typename Model::lattice_type;
  using RCluster = typename RClusterDmn::parameter_type;
  using KCluster = typename domains::cluster_symmetry<RCluster>::dual_type;
  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;
  using H0Type =
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, func::dmn_0<KCluster>>>;

  if constexpr (Lattice::DIMENSION != 2 ||
                !detail::LatticeHasInitializeH0<Lattice, Parameters, H0Type>::value) {
    (void)parameters;
    return "";
  }
  else {
    using KSymDmn = typename domains::cluster_symmetry<KCluster>::sym_super_cell_dmn_t;
    using Element = domains::point_group_symmetry_element;

    // The realized declared group, populated by the initializer that just ran. Only the first
    // get_size() entries are the live group: on a re-initialization (some tests run update_domains
    // twice per process) the search zeroes the count but appends to the never-cleared element
    // vector, so the vector can be longer than the group. The reference tracks the singleton
    // through the group installs below; the copy preserves the declared ops.
    const std::vector<Element>& realized = KSymDmn::parameter_type::get_elements();
    const int n_declared = KSymDmn::parameter_type::get_size();
    const std::vector<Element> declared(realized.begin(), realized.begin() + n_declared);

    // Derive: candidate pool -> same geometric filter -> H0-invariance gate.
    detail::installPointGroup<KCluster, domains::holohedry_pool_2D>();

    H0Type H0;
    Model::initializeH0(parameters, H0);

    std::vector<Element> derived;
    for (const int s : verifiedSymmetryOps<KCluster>(H0))
      derived.push_back(realized[s]);

    // Restore the declared production state, and verify the restoration: every consumer of the
    // symmetry domains must see exactly the state the declared initializer produced.
    detail::installPointGroup<KCluster, typename Lattice::DCA_point_group>();

    if (KSymDmn::parameter_type::get_size() != n_declared)
      throw std::logic_error(
          "deriveAndComparePointGroup: failed to restore the declared symmetry state of " +
          KCluster::get_name() + ".");
    for (std::size_t i = 0; i < declared.size(); ++i)
      if (!detail::isSameSymmetryElement(realized[i], declared[i]))
        throw std::logic_error(
            "deriveAndComparePointGroup: failed to restore the declared symmetry state of " +
            KCluster::get_name() + ".");

    // Compare as sets and report the divergence, if any.
    int n_underdeclared = 0;
    for (const auto& op : derived)
      if (!detail::containsSymmetryElement(declared, op))
        ++n_underdeclared;

    int n_unverified = 0;
    for (const auto& op : declared)
      if (!detail::containsSymmetryElement(derived, op))
        ++n_unverified;

    if (n_underdeclared == 0 && n_unverified == 0)
      return "";

    std::ostringstream report;
    report << "derived-symmetry check [" << KCluster::get_name() << "]:\n"
           << "\tdeclared group: " << declared.size() << " op(s); derived group: " << derived.size()
           << " op(s).\n";
    if (n_underdeclared > 0)
      report << "\tunder-declared: " << n_underdeclared
             << " derived H0-verified symmetry op(s) beyond the declared group "
             << "(symmetry averaging left unexploited).\n";
    if (n_unverified > 0)
      report << "\tWARNING: " << n_unverified
             << " declared op(s) not reproduced by the derivation -- the op fails the "
             << "H0-invariance check (or, the lattice is not in the pool's standard orientation).\n";
    return report.str();
  }
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_SYMMETRIZATION_DERIVE_POINT_GROUP_HPP
