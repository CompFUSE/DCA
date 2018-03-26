// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_HPP

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace kernelinterface{
// dca::phys::solver::ctint::kernelinterface::

double g0_interpolation_test(const double tau, const double lindex,
                     const double beta, const double n_div_beta,
                     const double* g0_minus, const double* values, const int* step_sbdm);

void compute_D(double* S, int lds, double* Q, int ldq, double* R, int ldr, int offset_i,
                      int offset_j, int n, const int* site_diff, int ldd, const int* sbdm_step,
                      const double* tau_r, const double* tau_l, const int* r_l, const int* r_r,
                      const int* nu_l, const int* nu_r, const ushort* aux_spin, int thread_id,
                      int stream_id);

}  // kernelinterface
}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_HPP
