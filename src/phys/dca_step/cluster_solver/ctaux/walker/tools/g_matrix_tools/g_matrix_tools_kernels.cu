// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements g_matrix_tools_kernels.hpp.

#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g_matrix_tools/g_matrix_tools_kernels.hpp"
#include <cassert>
#include "cuda_runtime.h"
#include "dca/linalg/util/error_cuda.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace gkernels {
// dca::phys::solver::ctaux::gkernels::

__global__ void read_G_matrix_kernel(int S, int vertex_index, int* i_index, int* j_index,
                                     bool* is_Bennett, double* exp_Vj, double* N_ptr, int LD_N,
                                     double* G_ptr, int LD_G, double* result_ptr, int incr) {
  int l = threadIdx.x;

  double result, delta;

  if (j_index[l] < vertex_index) {
    delta = i_index[l] == j_index[l] ? 1. : 0.;
    result = (N_ptr[i_index[l] + LD_N * j_index[l]] * exp_Vj[l] - delta) / (exp_Vj[l] - 1.);
  }
  else
    result = G_ptr[i_index[l] + LD_G * (j_index[l] - vertex_index)];

  result_ptr[l * incr] = is_Bennett[l] ? 0. : result;
}

void read_G_matrix_elements(int S, int vertex_index, int* i_ptr, int* j_ptr, bool* is_Bennett_ptr,
                            double* exp_Vj_ptr, double* N_ptr, int N_LD, double* G_ptr, int G_LD,
                            double* result_ptr, int incr) {
  read_G_matrix_kernel<<<1, S>>>(S, vertex_index, i_ptr, j_ptr, is_Bennett_ptr, exp_Vj_ptr, N_ptr,
                                 N_LD, G_ptr, G_LD, result_ptr, incr);

  checkErrorsCudaDebug();
}

__global__ void compute_row_on_Gamma_matrix_kernel(int row_index, int vertex_index, int* indices,
                                                   double* exp_V, double* N_ptr, int LD_N,
                                                   double* G_ptr, int LD_G, double* row_ptr,
                                                   int incr) {
  // int l = threadIdx.x;
  int l = blockIdx.x;

  int i_index, j_index;
  double delta;

  i_index = indices[row_index];
  j_index = indices[l];

  if (j_index < vertex_index) {
    delta = i_index == j_index ? 1 : 0;
    row_ptr[l * incr] = (N_ptr[i_index + LD_N * j_index] * exp_V[l] - delta) / (exp_V[l] - 1.);
  }
  else
    row_ptr[l * incr] = G_ptr[i_index + LD_G * (j_index - vertex_index)];
}

void compute_row_on_Gamma_matrix(int row_index, int S, int vertex_index, int* indices,
                                 double* exp_V, double* N_ptr, int LD_N, double* G_ptr, int LD_G,
                                 double* row_ptr, int incr) {
  // assert(S<get_number_of_threads());
  compute_row_on_Gamma_matrix_kernel<<<S, 1>>>(row_index, vertex_index, indices, exp_V, N_ptr, LD_N,
                                               G_ptr, LD_G, row_ptr, incr);

  checkErrorsCudaDebug();
}

__global__ void compute_col_on_Gamma_matrix_kernel(int col_index, int vertex_index, int* indices,
                                                   double* exp_V, double* N_ptr, int LD_N,
                                                   double* G_ptr, int LD_G, double* col_ptr,
                                                   int incr) {
  // int l = threadIdx.x;
  int l = blockIdx.x;

  int i_index, j_index;
  double delta, exp_Vj;

  i_index = indices[l];
  j_index = indices[col_index];

  exp_Vj = exp_V[col_index];

  if (j_index < vertex_index) {
    delta = i_index == j_index ? 1 : 0;
    col_ptr[l * incr] = (N_ptr[i_index + LD_N * j_index] * exp_Vj - delta) / (exp_Vj - 1.);
  }
  else
    col_ptr[l * incr] = G_ptr[i_index + LD_G * (j_index - vertex_index)];
}

void compute_col_on_Gamma_matrix(int col_index, int S, int vertex_index, int* indices,
                                 double* exp_V, double* N_ptr, int LD_N, double* G_ptr, int LD_G,
                                 double* col_ptr, int incr) {
  // assert(S<get_number_of_threads());
  compute_col_on_Gamma_matrix_kernel<<<S, 1>>>(col_index, vertex_index, indices, exp_V, N_ptr, LD_N,
                                               G_ptr, LD_G, col_ptr, incr);

  checkErrorsCudaDebug();
}

}  // gkernels
}  // ctaux
}  // solver
}  // phys
}  // dca
