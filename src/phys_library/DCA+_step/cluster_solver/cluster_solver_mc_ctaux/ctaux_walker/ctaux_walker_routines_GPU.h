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

#include "comp_library/linalg/src/matrix.h"
#include "comp_library/linalg/src/vector.h"

namespace DCA {
namespace QMCI {
namespace CT_AUX_WALKER_GPU_KERNELS {
// DCA::QMCI::CT_AUX_WALKER_GPU_KERNELS::

void compute_Gamma(double* Gamma, int Gamma_n, int Gamma_ld, double* N, int N_r, int N_c, int N_ld,
                   double* G, int G_r, int G_c, int G_ld, int* random_vertex_vector, double* exp_V,
                   double* exp_delta_V, int thread_id, int stream_id);

}  // CT_AUX_WALKER_GPU_KERNELS

template <>
class CT_AUX_WALKER_TOOLS<LIN_ALG::GPU> {
public:
  static void compute_Gamma(LIN_ALG::matrix<double, LIN_ALG::GPU>& Gamma,
                            LIN_ALG::matrix<double, LIN_ALG::GPU>& N,
                            LIN_ALG::matrix<double, LIN_ALG::GPU>& G,
                            LIN_ALG::vector<int, LIN_ALG::GPU>& random_vertex_vector,
                            LIN_ALG::vector<double, LIN_ALG::GPU>& exp_V,
                            LIN_ALG::vector<double, LIN_ALG::GPU>& exp_delta_V, int thread_id,
                            int stream_id) {
    Gamma.resize(random_vertex_vector.size());

    assert(Gamma.get_number_of_rows() == Gamma.get_number_of_cols());

    CT_AUX_WALKER_GPU_KERNELS::compute_Gamma(
        Gamma.get_ptr(), Gamma.get_number_of_rows(), Gamma.get_leading_dimension(), N.get_ptr(),
        N.get_number_of_rows(), N.get_number_of_cols(), N.get_leading_dimension(), G.get_ptr(),
        G.get_number_of_rows(), G.get_number_of_cols(), G.get_leading_dimension(),
        random_vertex_vector.get_ptr(), exp_V.get_ptr(), exp_delta_V.get_ptr(), thread_id, stream_id);
  }
};

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_ROUTINES_GPU_H
