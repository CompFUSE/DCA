// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Interface to the GPU kernels used by the TpAccumulator class.

#ifndef DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_KERNELS_INTERFACE_HPP
#define DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_KERNELS_INTERFACE_HPP

#include <complex>
#include <vector>

#include "dca/phys/four_point_type.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {
// dca::phys::solver::accumulator::details::

template <typename Real>
void computeGSingleband(std::complex<Real>* G, int ldg, const std::complex<Real>* G0, int nk,
                        int nw_pos, Real beta, cudaStream_t stream);

template <typename Real>
void computeGMultiband(std::complex<Real>* G, int ldg, const std::complex<Real>* G0, int ldg0,
                       int nb, int nk, int nw_pos, Real beta, cudaStream_t stream);

template <typename Real, FourPointType type>
void updateG4(std::complex<Real>* G4, const std::complex<Real>* G_up, int lggu,
              const std::complex<Real>* G_down, int ldgd, int nb, int nk, int nw_pos,
              int nw_exchange, int nk_exchange, int sign, cudaStream_t stream);

void initializeG4Helpers(int nb, int nk, int nw_pos, const std::vector<int>& delta_k_indices,
                         const std::vector<int>& delta_w_indices, const int* add_k, int lda,
                         const int* sub_k, int lds);

}  // details
}  // accumulator
}  // solver
}  // phys
}  // dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_KERNELS_INTERFACE_HPP
