// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// GPU kernels for walker tools.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_KERNELS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_KERNELS_HPP

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace walkerkernels {
// dca::phys::solver::ctaux::walkerkernels::

void compute_Gamma(double* Gamma, int Gamma_n, int Gamma_ld, double* N, int N_r, int N_c, int N_ld,
                   double* G, int G_r, int G_c, int G_ld, int* random_vertex_vector, double* exp_V,
                   double* exp_delta_V, int thread_id, int stream_id);

}  // walkerkernels
}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_CT_AUX_WALKER_TOOLS_KERNELS_HPP
