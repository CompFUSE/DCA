// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class manages the symmetry of the cluster.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_SYMMETRY_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_SYMMETRY_HPP

#include "dca/platform/dca_gpu.h"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_family.hpp"
#include "dca/phys/domains/cluster/dual_cluster.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/point_group_symmetry_domain.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
class cluster_symmetry<cluster_domain<scalar_type, D, N, R, S>> {
  const static int DIMENSION = D;

  const static CLUSTER_NAMES NAME = N;
  const static CLUSTER_REPRESENTATION REPRESENTATION = R;
  const static CLUSTER_SHAPE SHAPE = S;

  const static CLUSTER_REPRESENTATION DUAL_REPRESENTATION =
      dual_cluster<REPRESENTATION>::REPRESENTATION;

public:
  typedef cluster_domain_family<scalar_type, D, N, S> cluster_family_type;

  typedef cluster_domain<scalar_type, D, N, REPRESENTATION, S> this_type;
  typedef cluster_domain<scalar_type, D, N, DUAL_REPRESENTATION, S> dual_type;

  typedef func::dmn_0<domains::electron_band_domain> b_dmn_t;
  typedef func::dmn_0<this_type> c_dmn_t;

  typedef func::dmn_0<domains::point_group_symmetry_domain<domains::UNIT_CELL, cluster_family_type>>
      sym_unit_cell_dmn_t;
  typedef func::dmn_0<domains::point_group_symmetry_domain<domains::SUPER_CELL, cluster_family_type>>
      sym_super_cell_dmn_t;

  typedef func::dmn_variadic<func::dmn_variadic<c_dmn_t, b_dmn_t>, sym_super_cell_dmn_t> symmetry_matrix_dmn_t;

  // (band_row, band_col) x super-cell op: the orbital operation matrix U_S, one per op.
  typedef func::dmn_variadic<func::dmn_variadic<b_dmn_t, b_dmn_t>, sym_super_cell_dmn_t> orbital_op_dmn_t;

  // (cluster point) x super-cell op: the band-independent image of each point under each op -- today's
  // get_symmetry_matrix().first promoted out of the pair. (Momentum images do not depend on the band.)
  typedef func::dmn_variadic<c_dmn_t, sym_super_cell_dmn_t> mapped_point_dmn_t;

  static func::function<std::pair<int, int>, symmetry_matrix_dmn_t>& get_symmetry_matrix() {
    static func::function<std::pair<int, int>, symmetry_matrix_dmn_t> symmetry_matrix(
        "symmetry_matrix_super_cell");
    return symmetry_matrix;
  }

  // Image of each cluster point under each super-cell op: the ".first" of get_symmetry_matrix(),
  // split into its own band-independent accessor. Populated for the momentum cluster in
  // set_symmetry_matrices, where the k-image (from linear_transform of k) does not depend on the
  // band. The real-space point image IS band-dependent -- intra-cell orbital offsets move with the
  // site, so it cannot be separated this way -- and is not stored here. M1 shadow accessor; the
  // uniform symmetry relation reads G(k) = U_S G(mapped_point) U_S^dagger.
  static func::function<int, mapped_point_dmn_t>& get_mapped_point() {
    static func::function<int, mapped_point_dmn_t> mapped_point("mapped_point_super_cell");
    return mapped_point;
  }

  // Diagonal folding phase V(G) = diag(e^{i G.a_alpha}) per (cluster point, band, super-cell op).
  // Populated alongside get_symmetry_matrix() in set_symmetry_matrices. Runtime guard rejects (throws)
  // any phase that is not +/-1, which holds for every shipped (symmorphic, point-symmetry-only) model.
  // \todo Lift the storage to a complex phase to support a genuine non-+/-1 fold.
  static func::function<double, symmetry_matrix_dmn_t>& get_fold_phase() {
    static func::function<double, symmetry_matrix_dmn_t> fold_phase("fold_phase_super_cell");
    return fold_phase;
  }

  // Orbital operation U_S per super-cell op: the band x band matrix s.t. the symmetry relation is
  // G(k) = U_S G(S^-1 k) U_S^dagger. M1 shadow accessor; populated by a model-aware step (it needs
  // the Hamiltonian / the model's sign data, neither visible in this cluster-domain layer). v1 stores
  // it real -- the values are a signed permutation (+/-1) for every shipped symmorphic model; complex
  // U_S to support genuine orbital mixing is a future generalization.
  static func::function<double, orbital_op_dmn_t>& get_orbital_op() {
    static func::function<double, orbital_op_dmn_t> orbital_op("orbital_op_super_cell");
    return orbital_op;
  }
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_SYMMETRY_HPP
