// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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

template<class T>
void compute_G_cols(int N_i, int N_r, int N_c, const int* p_ptr, const T* exp_V_ptr, const T* N_ptr,
                    int N_ld, const T* G_ptr, int G_ld, T* G_cols_ptr, int G_cols_ld,
                    int thread_id, int stream_id);

template<class T>
void compute_d_vector(int N_i, const int* d_ind, T* d_ptr, int* p_ptr, const T* N_ptr, int N_ld,
                      int thread_id, int stream_id);
}  // nkernels
}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_N_MATRIX_TOOLS_N_MATRIX_TOOLS_KERNELS_HPP
