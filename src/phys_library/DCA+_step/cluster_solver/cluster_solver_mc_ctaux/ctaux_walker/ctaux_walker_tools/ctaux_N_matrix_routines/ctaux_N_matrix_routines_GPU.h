// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization for GPU.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_N_MATRIX_ROUTINES_CTAUX_N_MATRIX_ROUTINES_GPU_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_N_MATRIX_ROUTINES_CTAUX_N_MATRIX_ROUTINES_GPU_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_N_matrix_routines/ctaux_N_matrix_routines_TEM.h"

#include <cassert>
#include <vector>

#include "comp_library/linalg/linalg.hpp"

namespace DCA {
namespace QMCI {
namespace N_MATRIX_TOOLS_GPU_KERNELS {
// DCA::QMCI::N_MATRIX_TOOLS_GPU_KERNELS::

void compute_G_cols(int N_i, int N_r, int N_c, int* p_ptr, double* exp_V_ptr, double* N_ptr,
                    int N_ld, double* G_ptr, int G_ld, double* G_cols_ptr, int G_cols_ld,
                    int thread_id, int stream_id);

void compute_d_vector(int N_i, int* d_ind, double* d_ptr, int* p_ptr, double* N_ptr, int N_ld,
                      int thread_id, int stream_id);

}  // N_MATRIX_TOOLS_GPU_KERNELS

template <typename parameters_type>
class N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type> {
  const static int MAX_VERTEX_SINGLETS = 4;

  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_t;

public:
  N_MATRIX_TOOLS(int id, parameters_type& parameters_ref);

  int* get_permutation();
  void set_permutation(std::vector<int>& p);

  void set_d_vector(std::vector<int>& d_index, dca::linalg::Matrix<double, dca::linalg::GPU>& N,
                    dca::linalg::Vector<double, dca::linalg::CPU>& d_inv);

  void set_d_vector(dca::linalg::Vector<double, dca::linalg::CPU>& d_inv);

  void scale_rows(dca::linalg::Matrix<double, dca::linalg::GPU>& N);

  double* get_device_ptr(dca::linalg::Vector<double, dca::linalg::CPU>& v);

  void copy_rows(dca::linalg::Matrix<double, dca::linalg::GPU>& N,
                 dca::linalg::Matrix<double, dca::linalg::GPU>& N_new_spins);

  void compute_G_cols(std::vector<double>& exp_V, dca::linalg::Matrix<double, dca::linalg::GPU>& N,
                      dca::linalg::Matrix<double, dca::linalg::GPU>& G,
                      dca::linalg::Matrix<double, dca::linalg::GPU>& G_cols);

private:
  int thread_id;
  int stream_id;

  parameters_type& parameters;
  concurrency_type& concurrency;

  dca::linalg::Vector<int, dca::linalg::GPU> identity;
  dca::linalg::Vector<int, dca::linalg::GPU> permutation;

  dca::linalg::Vector<double, dca::linalg::GPU> tmp;

  dca::linalg::Vector<double, dca::linalg::GPU> exp_V;

  dca::linalg::Vector<int, dca::linalg::GPU> d_ind;
  dca::linalg::Vector<double, dca::linalg::GPU> d_vec;
};

template <typename parameters_type>
N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::N_MATRIX_TOOLS(int id,
                                                                  parameters_type& parameters_ref)
    : thread_id(id),
      stream_id(0),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      identity("identity    N_MATRIX_TOOLS<dca::linalg::GPU>",
               MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),
      permutation("permutation N_MATRIX_TOOLS<dca::linalg::GPU>",
                  MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),

      tmp("tmp   N_MATRIX_TOOLS<dca::linalg::GPU>", MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),
      exp_V("exp_V N_MATRIX_TOOLS<dca::linalg::GPU>", MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),

      d_ind("d_ind N_MATRIX_TOOLS<dca::linalg::GPU>", MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),
      d_vec("d_vec N_MATRIX_TOOLS<dca::linalg::GPU>", MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()) {
  {
    identity.setThreadAndStreamId(thread_id, stream_id);
    permutation.setThreadAndStreamId(thread_id, stream_id);

    tmp.setThreadAndStreamId(thread_id, stream_id);
    exp_V.setThreadAndStreamId(thread_id, stream_id);

    d_ind.setThreadAndStreamId(thread_id, stream_id);
    d_vec.setThreadAndStreamId(thread_id, stream_id);
  }

  {
    std::vector<int> id_tmp(MAX_VERTEX_SINGLETS * parameters.get_K_PHANI());

    for (int l = 0; l < MAX_VERTEX_SINGLETS * parameters.get_K_PHANI(); ++l)
      id_tmp[l] = l;

    identity.set(id_tmp);
  }
}

template <typename parameters_type>
int* N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::get_permutation() {
  return permutation.ptr();
}

template <typename parameters_type>
void N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::set_permutation(std::vector<int>& p) {
  permutation.set(p, LIN_ALG::ASYNCHRONOUS);
}

template <typename parameters_type>
void N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::set_d_vector(
    dca::linalg::Vector<double, dca::linalg::CPU>& d_inv) {
  d_vec.set(d_inv, LIN_ALG::ASYNCHRONOUS);
}

template <typename parameters_type>
void N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::scale_rows(
    dca::linalg::Matrix<double, dca::linalg::GPU>& N) {
  assert(permutation.size() == d_vec.size());

  int N_i = permutation.size();
  int N_c = N.nrCols();

  int N_LD = N.leadingDimension();

  LIN_ALG::SCALE<dca::linalg::GPU>::many_rows(N_c, N_i, permutation.ptr(), d_vec.ptr(), N.ptr(),
                                              N_LD, thread_id, stream_id);
}

template <typename parameters_type>
double* N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::get_device_ptr(
    dca::linalg::Vector<double, dca::linalg::CPU>& v) {
  tmp.set(v, LIN_ALG::ASYNCHRONOUS);

  return tmp.ptr();
}

template <typename parameters_type>
void N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::copy_rows(
    dca::linalg::Matrix<double, dca::linalg::GPU>& N,
    dca::linalg::Matrix<double, dca::linalg::GPU>& N_new_spins) {
  assert(N_new_spins.nrCols() == N.nrCols());
  assert(N_new_spins.nrRows() == permutation.size());

  int N_i = permutation.size();
  int N_c = N.nrCols();

  assert(N_i <= identity.size());

  LIN_ALG::COPY<dca::linalg::GPU>::many_rows(N_c, N_i, permutation.ptr(), N.ptr(),
                                             N.leadingDimension(), identity.ptr(), N_new_spins.ptr(),
                                             N_new_spins.leadingDimension(), thread_id, stream_id);
}

template <typename parameters_type>
void N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::compute_G_cols(
    std::vector<double>& exp_V_CPU, dca::linalg::Matrix<double, dca::linalg::GPU>& N,
    dca::linalg::Matrix<double, dca::linalg::GPU>& G,
    dca::linalg::Matrix<double, dca::linalg::GPU>& G_cols) {
  exp_V.set(exp_V_CPU, LIN_ALG::ASYNCHRONOUS);

  assert(N.nrRows() == G.nrRows());
  assert(N.nrRows() == G_cols.nrRows());
  assert(exp_V.size() == permutation.size());

  int N_i = permutation.size();
  int N_r = N.nrRows();

  int N_c = N.nrCols() - G.nrCols();

  N_MATRIX_TOOLS_GPU_KERNELS::compute_G_cols(
      N_i, N_r, N_c, permutation.ptr(), exp_V.ptr(), N.ptr(), N.leadingDimension(), G.ptr(),
      G.leadingDimension(), G_cols.ptr(), G_cols.leadingDimension(), thread_id, stream_id);
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_N_MATRIX_ROUTINES_CTAUX_N_MATRIX_ROUTINES_GPU_H
