// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// GPU kernels for walker tools.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_KERNELS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_KERNELS_HPP

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace walkerkernels {
// dca::phys::solver::ctaux::walkerkernels::

template <typename Real>
void compute_Gamma(Real* Gamma, int Gamma_n, int Gamma_ld, Real* N, int N_r, int N_c, int N_ld,
                   Real* G, int G_r, int G_c, int G_ld, int* random_vertex_vector, Real* exp_V,
                   Real* exp_delta_V, int thread_id, int stream_id);

}  // namespace walkerkernels
}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_KERNELS_HPP
