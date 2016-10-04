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

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G_MATRIX_ROUTINES_CTAUX_G_MATRIX_ROUTINES_GPU_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G_MATRIX_ROUTINES_CTAUX_G_MATRIX_ROUTINES_GPU_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_G_matrix_routines/ctaux_G_matrix_routines_TEM.h"

#include <cassert>

#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

namespace DCA {
namespace QMCI {
namespace GPU_KERNELS_G_TOOLS {
// DCA::QMCI::GPU_KERNELS_G_TOOLS::

void read_G_matrix_elements(int S, int vertex_index, int* i_index, int* j_index, bool* is_Bennett,
                            double* exp_Vj, double* N, int LD_N, double* G_precomputed, int LD_G,
                            double* result_ptr, int incr);

void compute_row_on_Gamma_matrix(int row_index, int S, int vertex_index, int* indices,
                                 double* exp_V, double* N, int LD_N, double* G_precomputed,
                                 int LD_G, double* row_ptr, int incr);

void compute_col_on_Gamma_matrix(int row_index, int S, int vertex_index, int* indices,
                                 double* exp_V, double* N, int LD_N, double* G_precomputed,
                                 int LD_G, double* row_ptr, int incr);

}  // GPU_KERNELS_G_TOOLS

template <typename parameters_type>
class G_MATRIX_TOOLS<dca::linalg::GPU, parameters_type> {
  const static int MAX_VERTEX_SINGLETS = 4;

  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_t;

public:
  G_MATRIX_TOOLS(int id, parameters_type& parameters_ref)
      : thread_id(id),

        parameters(parameters_ref),
        concurrency(parameters.get_concurrency()),

        i_index_gpu(0, MAX_VERTEX_SINGLETS * parameters_ref.get_submatrix_size()),
        j_index_gpu(0, MAX_VERTEX_SINGLETS * parameters_ref.get_submatrix_size()),
        is_Bennett_gpu(0, MAX_VERTEX_SINGLETS * parameters_ref.get_submatrix_size()),
        exp_Vj_gpu(0, MAX_VERTEX_SINGLETS * parameters_ref.get_submatrix_size()) {}

  void read_G_matrix_elements(dca::linalg::Vector<int, dca::linalg::CPU>& i_index,
                              dca::linalg::Vector<int, dca::linalg::CPU>& j_index,
                              dca::linalg::Vector<bool, dca::linalg::CPU>& is_Bennett,
                              dca::linalg::Vector<double, dca::linalg::CPU>& exp_Vj,
                              dca::linalg::Matrix<double, dca::linalg::GPU>& N,
                              dca::linalg::Matrix<double, dca::linalg::GPU>& G_precomputed,
                              double* result_ptr, int incr) {
    assert(i_index.size() == j_index.size());
    assert(i_index.size() == is_Bennett.size());
    assert(i_index.size() == exp_Vj.size());

    i_index_gpu = i_index;
    j_index_gpu = j_index;
    is_Bennett_gpu = is_Bennett;
    exp_Vj_gpu = exp_Vj;

    int vertex_index = N.nrCols() - G_precomputed.nrCols();

    GPU_KERNELS_G_TOOLS::read_G_matrix_elements(
        i_index_gpu.size(), vertex_index, i_index_gpu.ptr(), j_index_gpu.ptr(),
        is_Bennett_gpu.ptr(), exp_Vj_gpu.ptr(), N.ptr(), N.leadingDimension(), G_precomputed.ptr(),
        G_precomputed.leadingDimension(), result_ptr, incr);
  }

  void compute_row_on_Gamma_matrix(int row_index, dca::linalg::Vector<int, dca::linalg::GPU>& indices,
                                   dca::linalg::Vector<double, dca::linalg::GPU>& exp_V,
                                   dca::linalg::Matrix<double, dca::linalg::GPU>& N,
                                   dca::linalg::Matrix<double, dca::linalg::GPU>& G_precomputed,
                                   double* row_ptr, int incr) {
    int vertex_index = N.nrCols() - G_precomputed.nrCols();

    GPU_KERNELS_G_TOOLS::compute_row_on_Gamma_matrix(
        row_index, indices.size(), vertex_index,

        indices.ptr(), exp_V.ptr(),

        N.ptr(), N.leadingDimension(), G_precomputed.ptr(), G_precomputed.leadingDimension(),
        row_ptr, incr);
  }

  void compute_col_on_Gamma_matrix(int col_index, dca::linalg::Vector<int, dca::linalg::GPU>& indices,
                                   dca::linalg::Vector<double, dca::linalg::GPU>& exp_V,
                                   dca::linalg::Matrix<double, dca::linalg::GPU>& N,
                                   dca::linalg::Matrix<double, dca::linalg::GPU>& G_precomputed,
                                   double* col_ptr, int incr) {
    int vertex_index = N.nrCols() - G_precomputed.nrCols();

    GPU_KERNELS_G_TOOLS::compute_col_on_Gamma_matrix(
        col_index, indices.size(), vertex_index,

        indices.ptr(), exp_V.ptr(),

        N.ptr(), N.leadingDimension(), G_precomputed.ptr(), G_precomputed.leadingDimension(),
        col_ptr, incr);
  }

private:
  int thread_id;

  parameters_type& parameters;
  concurrency_type& concurrency;

  dca::linalg::Vector<int, dca::linalg::GPU> i_index_gpu, j_index_gpu;
  dca::linalg::Vector<bool, dca::linalg::GPU> is_Bennett_gpu;
  dca::linalg::Vector<double, dca::linalg::GPU> exp_Vj_gpu;
};

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G_MATRIX_ROUTINES_CTAUX_G_MATRIX_ROUTINES_GPU_H
