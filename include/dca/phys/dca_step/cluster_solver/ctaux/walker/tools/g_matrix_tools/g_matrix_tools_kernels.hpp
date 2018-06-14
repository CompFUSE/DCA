// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// GPU kernels for G-matrix tools.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G_MATRIX_TOOLS_G_MATRIX_TOOLS_KERNELS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G_MATRIX_TOOLS_G_MATRIX_TOOLS_KERNELS_HPP

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace gkernels {
// dca::phys::solver::ctaux::gkernels::

void read_G_matrix_elements(int S, int vertex_index, int* i_ptr, int* j_ptr, bool* is_Bennett_ptr,
                            double* exp_Vj_ptr, double* N_ptr, int N_LD, double* G_ptr, int G_LD,
                            double* result_ptr, int incr);

void compute_row_on_Gamma_matrix(int row_index, int S, int vertex_index, int* indices,
                                 double* exp_V, double* N_ptr, int LD_N, double* G_ptr, int LD_G,
                                 double* row_ptr, int incr);

void compute_col_on_Gamma_matrix(int col_index, int S, int vertex_index, int* indices,
                                 double* exp_V, double* N_ptr, int LD_N, double* G_ptr, int LD_G,
                                 double* col_ptr, int incr);

}  // gkernels
}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G_MATRIX_TOOLS_G_MATRIX_TOOLS_KERNELS_HPP
