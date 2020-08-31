// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// GPU kernels for G0 interpolation.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_KERNELS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_KERNELS_HPP

#include <utility>

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace g0kernels {
// dca::phys::solver::ctaux::g0kernels::

template <typename Real>
void interpolate_G0_matrix_on_GPU(int Nb, int Nr, int Nt, double beta, int Nv, int* b, int* r,
                                  Real* t, Real* G0, std::pair<int, int> G0_cs,
                                  std::pair<int, int> G0_gs, Real* r0_min_r1,
                                  std::pair<int, int> r0_min_r1_cs,
                                  std::pair<int, int> r0_min_r1_gs, Real* G0_r_t,
                                  std::pair<int, int> G0_r_t_cs, std::pair<int, int> G0_r_t_gs,
                                  Real* grad_G0_r_t, std::pair<int, int> grad_G0_r_t_cs,
                                  std::pair<int, int> grad_G0_r_t_gs);

template <typename Real>
void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, Real beta, int Nc, int Nv, int* b, int* r,
                                Real* t, Real* G0, std::pair<int, int> G0_cs,
                                std::pair<int, int> G0_gs, Real* r0_min_r1,
                                std::pair<int, int> r0_min_r1_cs, std::pair<int, int> r0_min_r1_gs,
                                Real* alpha, std::pair<int, int> alpha_cs,
                                std::pair<int, int> alpha_gs);

template <typename Real>
void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, Real beta, int Nc, int Nv, int* b, int* r,
                                Real* t, Real* G0, std::pair<int, int> G0_cs,
                                std::pair<int, int> G0_gs, Real* r0_min_r1,
                                std::pair<int, int> r0_min_r1_cs, std::pair<int, int> r0_min_r1_gs,
                                Real* alpha, std::pair<int, int> alpha_cs,
                                std::pair<int, int> alpha_gs, int thread_id, int stream_id);

}  // namespace g0kernels
}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_KERNELS_HPP
