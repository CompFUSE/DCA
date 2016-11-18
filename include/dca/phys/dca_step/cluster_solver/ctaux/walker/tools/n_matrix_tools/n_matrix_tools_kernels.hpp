// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// GPU kernels for N-matrix tools.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_N_MATRIX_TOOLS_N_MATRIX_TOOLS_KERNELS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_N_MATRIX_TOOLS_N_MATRIX_TOOLS_KERNELS_HPP

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace nkernels {
// dca::phys::solver::ctaux::nkernels::

void compute_G_cols(int N_i, int N_r, int N_c, int* p_ptr, double* exp_V_ptr, double* N_ptr,
                    int N_ld, double* G_ptr, int G_ld, double* G_cols_ptr, int G_cols_ld,
                    int thread_id, int stream_id);

void compute_d_vector(int N_i, int* d_ind, double* d_ptr, int* p_ptr, double* N_ptr, int N_ld,
                      int thread_id, int stream_id);
}  // nkernels
}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_N_MATRIX_TOOLS_N_MATRIX_TOOLS_KERNELS_HPP
