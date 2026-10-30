// Copyright (C) 2026 UT-Battelle, LLC
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
// Scope: for every analytic model, U_S is a real signed permutation U_S = D P. Both halves are
// derived from H0: the permutation P (which orbital maps to which) and the diagonal of +/-1 signs
// D = diag(sigma). The symmetry table's band image (.second) is only the first candidate for P;
// when H0 rejects it -- geometry cannot see an on-site orbital rotation, such as the dxz/dyz swap
// under C4 of two orbitals sharing a_vec = 0 -- the remaining permutations are searched. n_b <= 3
// for every shipped model, so that is at most 3! = 6 candidates.
//
// The symmetry table stores not the exact momentum image S k but its folded representative
// k_new = S k - G, and with intra-cell orbital offsets a_b the fold is not free: H0 picks up the
// diagonal folding phase V(G) = diag(phi_b), phi_b = e^{i G.a_b} (cluster_symmetry::get_fold_phase(),
// a real +/-1 for every shipped model). Substituting U_S = D P into the symmetry relation
// H0(S k) = U_S H0(k) U_S^dagger and folding S k back to k_new gives, entry by entry,
//
//   sigma(b0) sigma(b1) H0(k)(b0,b1) = phi(image(b0)) phi(image(b1)) H0(k_new)(image(b0), image(b1)).
//
// So each nonzero coupling fixes the product sigma(b0) sigma(b1) to the +/-1 fold-corrected ratio
// of the two H0 entries. The phi factors matter: without them the solve absorbs the folding phase
// into the signs and returns a representative-convention-dependent dressing of U_S (e.g. masking
// the odd parity of p_x under a mirror at a zone-boundary momentum) instead of the intrinsic
// orbital transform.
//
// For a candidate permutation, the signs are found by enumerating the at most 2^(n_b-1) <= 4
// vectors with the gauge sigma(0) = +1 and keeping the one that satisfies every constraint.
// A malformed candidate (an out-of-range or repeated band image -- only the geometric one can be)
// is rejected up front; an admissible one is rejected, and the next tried, when:
//
//    * a coupling vanishes on only one side (a magnitude mismatch -- the permutation is not a
//      symmetry of H0)
//    * a fold-corrected ratio is not a real +/-1 (this permutation would need a
//      non-signed-permutation U_S)
//    * no sign vector satisfies all couplings
//
// Only when every band permutation is rejected does the operation fail outright: it then needs a
// genuinely non-signed-permutation U_S (orbital mixing or a complex phase out of scope) or is not
// a symmetry of H0.

#ifndef DCA_PHYS_DCA_STEP_SYMMETRIZATION_SOLVE_ORBITAL_OP_SIGNS_HPP
#define DCA_PHYS_DCA_STEP_SYMMETRIZATION_SOLVE_ORBITAL_OP_SIGNS_HPP

#include <algorithm>
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

// Tests one candidate band permutation `image` against H0. Returns true, and writes that op into
// u_s, iff `image` is a signed-permutation symmetry of H0: every nonzero coupling gives a real
// +/-1 (fold-corrected) ratio, and some sign vector satisfies the resulting products.
template <typename SymFunc, typename FoldFunc, typename H0Function, typename UFunc>
bool tryPermutation(const std::vector<int>& image, int s, int nb, int nk, const SymFunc& sym,
                    const FoldFunc& fold_phase, const H0Function& H0, UFunc& u_s) {
  // Reject a malformed candidate before it indexes H0: an out-of-range entry (the -1 recorded when
  // set_symmetry_matrices finds no admissible image) or a repeated one (two orbitals sharing a site
  // and a flavor collide).
  std::vector<bool> is_hit(nb, false);
  for (int b = 0; b < nb; ++b) {
    if (image[b] < 0 || image[b] >= nb || is_hit[image[b]])
      return false;
    is_hit[image[b]] = true;
  }

  // Gather one sign-product constraint per nonzero coupling, over all cluster momenta. Each nonzero
  // coupling forces sigma(b0) sigma(b1) to the +/-1 fold-corrected ratio of the two H0 entries; the
  // ratio must be a real +/-1 for a signed permutation to exist under this permutation at all.
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
          return false;  // magnitude mismatch: this permutation is not a symmetry of H0.
        // Undress the fold before taking the ratio. phi lives at the image bands -- the rhs entry
        // sits at (image(b0), image(b1)) -- with G determined by (k, s); the stored table carries
        // the band as a free index, so the image-band lookup is direct.
        const double fold = fold_phase(k, image[b0], s) * fold_phase(k, image[b1], s);
        const std::complex<double> ratio = lhs / (fold * rhs);
        if (std::abs(std::imag(ratio)) > kSignRealTol ||
            std::abs(std::abs(ratio) - 1.) > kSignRealTol)
          return false;  // non-signed-permutation ratio: this permutation would need a dense U_S.

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
  if (!found)
    return false;  // no sign vector reproduces H0 under this permutation.

  // Materialize U_S = D P one entry at a time: sign sigma(b) at (row b, column image[b]).
  for (int b = 0; b < nb; ++b)
    u_s(b, image[b], s) = static_cast<double>(sigma[b]);
  return true;
}

// Derives op s's orbital operation U_S = D P from H0 -- both the permutation P and the signs D --
// and writes it into u_s. H0 is the authority: U_S is whatever signed permutation makes
// H0(S k) = U_S H0(k) U_S^dagger hold.
template <typename SymFunc, typename FoldFunc, typename H0Function, typename UFunc>
void deriveOrbitalOpForOp(int s, int nb, int nk, const SymFunc& sym, const FoldFunc& fold_phase,
                          const H0Function& H0, UFunc& u_s) {
  // The geometric image is the canonical candidate: verify it against H0 first.
  std::vector<int> geometric(nb);
  for (int b = 0; b < nb; ++b)
    geometric[b] = sym(0, b, s).second;

  if (tryPermutation(geometric, s, nb, nk, sym, fold_phase, H0, u_s))
    return;

  // H0 did not admit it; check the remaining permutations, in std::next_permutation order so the
  // search is deterministic, skipping the geometric candidate already tried.
  std::vector<int> image(nb);
  for (int b = 0; b < nb; ++b)
    image[b] = b;
  do {
    if (image != geometric && tryPermutation(image, s, nb, nk, sym, fold_phase, H0, u_s))
      return;
  } while (std::next_permutation(image.begin(), image.end()));

  throw std::logic_error(
      "solveOrbitalOpSignsFromH0: no signed-permutation U_S reproduces H0 under op " +
      std::to_string(s) +
      " for any band permutation -- the operation needs a non-signed-permutation U_S (genuine "
      "orbital mixing or a complex phase) or is not a symmetry of H0.");
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
  const auto& fold_phase = domains::cluster_symmetry<KCluster>::get_fold_phase();
  auto& u_s = domains::cluster_symmetry<KCluster>::get_orbital_op();
  u_s = 0.;

  const int nb = BDmn::dmn_size();
  const int nk = CDmn::dmn_size();
  const int n_ops = SymDmn::dmn_size();

  for (int s = 0; s < n_ops; ++s)
    detail::deriveOrbitalOpForOp(s, nb, nk, sym, fold_phase, H0, u_s);
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
  const auto& fold_phase = domains::cluster_symmetry<KCluster>::get_fold_phase();
  auto& u_s = domains::cluster_symmetry<KCluster>::get_orbital_op();
  u_s = 0.;

  const int nb = BDmn::dmn_size();
  const int nk = CDmn::dmn_size();
  const int n_ops = SymDmn::dmn_size();

  std::vector<int> verified;
  for (int s = 0; s < n_ops; ++s) {
    try {
      detail::deriveOrbitalOpForOp(s, nb, nk, sym, fold_phase, H0, u_s);
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
