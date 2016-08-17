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

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_ROUTINES_GPU_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_ROUTINES_GPU_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_routines_template.h"

#include <cassert>

#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

namespace DCA {
namespace QMCI {
namespace CT_AUX_WALKER_GPU_KERNELS {
// DCA::QMCI::CT_AUX_WALKER_GPU_KERNELS::

void compute_Gamma(double* Gamma, int Gamma_n, int Gamma_ld, double* N, int N_r, int N_c, int N_ld,
                   double* G, int G_r, int G_c, int G_ld, int* random_vertex_vector, double* exp_V,
                   double* exp_delta_V, int thread_id, int stream_id);

}  // CT_AUX_WALKER_GPU_KERNELS

template <>
class CT_AUX_WALKER_TOOLS<dca::linalg::GPU> {
public:
  static void compute_Gamma(dca::linalg::Matrix<double, dca::linalg::GPU>& Gamma,
                            dca::linalg::Matrix<double, dca::linalg::GPU>& N,
                            dca::linalg::Matrix<double, dca::linalg::GPU>& G,
                            dca::linalg::Vector<int, dca::linalg::GPU>& random_vertex_vector,
                            dca::linalg::Vector<double, dca::linalg::GPU>& exp_V,
                            dca::linalg::Vector<double, dca::linalg::GPU>& exp_delta_V,
                            int thread_id, int stream_id) {
    Gamma.resize(random_vertex_vector.size());

    assert(Gamma.nrRows() == Gamma.nrCols());

    CT_AUX_WALKER_GPU_KERNELS::compute_Gamma(
        Gamma.ptr(), Gamma.nrRows(), Gamma.leadingDimension(), N.ptr(), N.nrRows(), N.nrCols(),
        N.leadingDimension(), G.ptr(), G.nrRows(), G.nrCols(), G.leadingDimension(),
        random_vertex_vector.ptr(), exp_V.ptr(), exp_delta_V.ptr(), thread_id, stream_id);
  }
};

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_ROUTINES_GPU_H
