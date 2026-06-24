// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Tyler Sax (tylersax@gmail.com)
//
// Derives the orbital-operation matrix U_S of each point-group operation from H0(k) and populates
// cluster_symmetry::get_orbital_op(), without reading hand-coded Lattice::transformationSignOf*.
//
// Scope: for every analytic model U_S is a real signed permutation U_S = D P, where the permutation
// P is already encoded in the symmetry table (the band image .second) and the only unknown is the
// diagonal of +/-1 signs D = diag(sigma). Substituting into the relation H0(k) = U_S H0(k_new) U_S^dagger
// gives, entry by entry,
//
//   H0(k)(b0,b1) = sigma(b0) sigma(b1) H0(k_new)(image(b0), image(b1)).
//
// So each nonzero coupling fixes the product sigma(b0) sigma(b1) to the +/-1 ratio of the two H0
// entries. The analytic models have at most n_b = 3 orbitals, so once the gauge sigma(0) = +1 is
// fixed there are at most 2^(n_b-1) <= 4 candidate sign vectors. We simply enumerate them and keep
// the one that satisfies every coupling constraint. Three conditions make an operation fail:
//
//   * a magnitude mismatch (an entry vanishing on only one side) means the permutation is not a
//     symmetry of H0;
//   * a non-real or non-unit ratio means the op needs a non-signed-permutation U_S (genuine orbital
//     mixing, not yet supported) or is not a symmetry;
//   * no candidate sign vector satisfying all couplings means no signed-permutation U_S exists.

#ifndef DCA_PHYS_DCA_STEP_SYMMETRIZATION_SOLVE_ORBITAL_OP_SIGNS_HPP
#define DCA_PHYS_DCA_STEP_SYMMETRIZATION_SOLVE_ORBITAL_OP_SIGNS_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_symmetry.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"

namespace dca {
namespace phys {
// dca::phys::

namespace detail {
// dca::phys::detail::

// Absolute threshold below which an H0 entry counts as a zero (no coupling).
constexpr double kCouplingZeroTol = 1e-10;
// Tolerance on the H0 ratio being a real +/-1: the imaginary part and the deviation of the
// magnitude from 1 must both be below this.
constexpr double kSignRealTol = 1e-9;

// Solves the +/-1 signs of one operation s and writes that op's U_S block into u_s. Throws if op s
// is not a signed-permutation symmetry of H0. Shared by the whole-group populator -- which lets the
// throw propagate, since a model must be a genuine symmetry group -- and by the group verification,
// which catches it to classify the op.
template <typename SymFunc, typename H0Function, typename UFunc>
void solveSignsForOp(int s, int nb, int nk, const SymFunc& sym, const H0Function& H0, UFunc& u_s) {
  // Band permutation P: row b's single entry lands in column image[b]. This is the geometry-derived
  // half of U_S (which orbital maps to which); it is k-independent for a point-group op, so read it
  // at k = 0. Only the signs that decorate it are unknown.
  std::vector<int> image(nb);
  for (int b = 0; b < nb; ++b)
    image[b] = sym(0, b, s).second;

  // Gather one sign-product constraint per nonzero coupling, over all cluster momenta. Each nonzero
  // coupling forces sigma(b0) sigma(b1) to the +/-1 ratio of the two H0 entries. The ratio must be a
  // real +/-1 for a signed permutation to exist at all; the two checks below (magnitude match, real
  // +/-1) are independent of the signs themselves, so they are tested here and reported as the first
  // two rejection branches.
  struct Constraint {
    int b0, b1, product;
  };
  std::vector<Constraint> constraints;
  for (int k = 0; k < nk; ++k) {
    const int k_new = sym(k, 0, s).first;
    for (int b0 = 0; b0 < nb; ++b0)
      for (int b1 = 0; b1 < nb; ++b1) {
        const std::complex<double> lhs = H0(b0, 0, b1, 0, k);
        const std::complex<double> rhs = H0(image[b0], 0, image[b1], 0, k_new);

        const bool zero_lhs = std::abs(lhs) < kCouplingZeroTol;
        const bool zero_rhs = std::abs(rhs) < kCouplingZeroTol;
        if (zero_lhs && zero_rhs)
          continue;  // coupling absent on both sides: consistent, but no constraint.
        if (zero_lhs != zero_rhs)
          throw std::logic_error(
              "solveOrbitalOpSignsFromH0: H0 magnitude mismatch for op " + std::to_string(s) +
              " -- the band permutation is not a symmetry of H0 (a coupling vanishes on only one "
              "side).");
        const std::complex<double> ratio = lhs / rhs;
        if (std::abs(std::imag(ratio)) > kSignRealTol ||
            std::abs(std::abs(ratio) - 1.) > kSignRealTol)
          throw std::domain_error(
              "solveOrbitalOpSignsFromH0: H0 ratio is not a real +/-1 for op " + std::to_string(s) +
              " -- the operation requires a non-signed-permutation U_S (genuine orbital mixing) or "
              "is not a symmetry of H0.");

        constraints.push_back({b0, b1, std::real(ratio) > 0. ? 1 : -1});
      }
  }

  // Find the signs by trying every assignment with the gauge sigma(0) = +1 -- 2^(n_b-1) <= 4
  // candidates for the supported models. A candidate is U_S iff it satisfies every gathered product
  // constraint; for a connected coupling graph (every analytic model) it is unique. The undetermined
  // global sign of any orbital that never couples is left at +1 by starting the search from the
  // all-plus assignment, matching the lowest-index-+1 gauge.
  std::vector<int> sigma(nb, 1);
  bool found = false;
  for (int mask = 0; mask < (1 << (nb - 1)) && !found; ++mask) {
    for (int b = 1; b < nb; ++b)
      sigma[b] = (mask & (1 << (b - 1))) ? -1 : 1;  // sigma[0] stays +1, the gauge.

    found = true;
    for (const Constraint& c : constraints)
      if (sigma[c.b0] * sigma[c.b1] != c.product) {
        found = false;
        break;
      }
  }
  // No sign vector satisfied every constraint, so no signed permutation reproduces H0 under this op
  if (!found)
    throw std::logic_error(
        "solveOrbitalOpSignsFromH0: no consistent +/-1 signs for op " + std::to_string(s) +
        " -- the H0 couplings admit no signed-permutation U_S.");

  // Materialize U_S = D P one entry at a time: sign sigma(b) at (row b, column image[b]).
  for (int b = 0; b < nb; ++b)
    u_s(b, image[b], s) = static_cast<double>(sigma[b]);
}

}  // namespace detail
// dca::phys::

// Populates cluster_symmetry<KCluster>::get_orbital_op() by solving each op's signs from H0.
// Throws if any op is not a signed-permutation symmetry of H0 (a model handed to this must be a
// genuine symmetry group).
template <typename KCluster, typename H0Function>
void solveOrbitalOpSignsFromH0(const H0Function& H0) {
  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SymDmn = typename domains::cluster_symmetry<KCluster>::sym_super_cell_dmn_t;
  using CDmn = typename domains::cluster_symmetry<KCluster>::c_dmn_t;

  const auto& sym = domains::cluster_symmetry<KCluster>::get_symmetry_matrix();
  auto& u_s = domains::cluster_symmetry<KCluster>::get_orbital_op();
  u_s = 0.;

  const int nb = BDmn::dmn_size();
  const int nk = CDmn::dmn_size();
  const int n_ops = SymDmn::dmn_size();

  for (int s = 0; s < n_ops; ++s)
    detail::solveSignsForOp(s, nb, nk, sym, H0, u_s);
}

// Returns the sorted op indices that pass the sign-consistency (H0-invariance) check -- the
// H0-verified symmetry group. Non-throwing: an op that fails any branch is simply omitted. The
// candidate ops come from the symmetry table, which today is the declared point group filtered by
// lattice geometry.
template <typename KCluster, typename H0Function>
std::vector<int> verifiedSymmetryOps(const H0Function& H0) {
  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SymDmn = typename domains::cluster_symmetry<KCluster>::sym_super_cell_dmn_t;
  using CDmn = typename domains::cluster_symmetry<KCluster>::c_dmn_t;

  const auto& sym = domains::cluster_symmetry<KCluster>::get_symmetry_matrix();
  auto& u_s = domains::cluster_symmetry<KCluster>::get_orbital_op();
  u_s = 0.;

  const int nb = BDmn::dmn_size();
  const int nk = CDmn::dmn_size();
  const int n_ops = SymDmn::dmn_size();

  std::vector<int> verified;
  for (int s = 0; s < n_ops; ++s) {
    try {
      detail::solveSignsForOp(s, nb, nk, sym, H0, u_s);
      verified.push_back(s);
    }
    catch (const std::exception&) {
      // op s is not a signed-permutation symmetry of H0 -> not part of the verified group.
    }
  }
  return verified;
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_SYMMETRIZATION_SOLVE_ORBITAL_OP_SIGNS_HPP
