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

void interpolate_G0_matrix_on_GPU(int Nb, int Nr, int Nt, double beta, int Nv, int* b, int* r,
                                  double* t, double* G0, std::pair<int, int> G0_cs,
                                  std::pair<int, int> G0_gs, double* r0_min_r1,
                                  std::pair<int, int> r0_min_r1_cs,
                                  std::pair<int, int> r0_min_r1_gs, double* G0_r_t,
                                  std::pair<int, int> G0_r_t_cs, std::pair<int, int> G0_r_t_gs,
                                  double* grad_G0_r_t, std::pair<int, int> grad_G0_r_t_cs,
                                  std::pair<int, int> grad_G0_r_t_gs);

void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, double beta, int Nc, int Nv, int* b, int* r,
                                double* t, double* G0, std::pair<int, int> G0_cs,
                                std::pair<int, int> G0_gs, double* r0_min_r1,
                                std::pair<int, int> r0_min_r1_cs, std::pair<int, int> r0_min_r1_gs,
                                double* alpha, std::pair<int, int> alpha_cs,
                                std::pair<int, int> alpha_gs);

void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, double beta, int Nc, int Nv, int* b, int* r,
                                double* t, double* G0, std::pair<int, int> G0_cs,
                                std::pair<int, int> G0_gs, double* r0_min_r1,
                                std::pair<int, int> r0_min_r1_cs, std::pair<int, int> r0_min_r1_gs,
                                double* alpha, std::pair<int, int> alpha_cs,
                                std::pair<int, int> alpha_gs, int thread_id, int stream_id);

}  // g0kernels
}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_KERNELS_HPP
