// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Interface to the GPU kernels used by cached_ndft_gpu.hpp.

#ifndef DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_KERNELS_INTERFACE_HPP
#define DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_KERNELS_INTERFACE_HPP

#include <complex>

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/triple.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {
// dca::phys::solver::accumulator::details::

template <typename InpScalar, typename Real>
void sortM(int size, const InpScalar* M, int ldm, std::complex<Real>* sorted_M, int lds,
           const Triple<Real>* config1, const Triple<Real>* config2, cudaStream_t stream);

template <typename Real>
void computeT(int n, int m, std::complex<Real>* T, int ldt, const Triple<Real>* config,
              const Real* w, bool transposed, cudaStream_t stream);

template <typename Real>
void rearrangeOutput(int nw, int no, int nb, const std::complex<Real>* in, int ldi,
                     std::complex<Real>* out, int ldo, cudaStream_t stream);

}  // namespace details
}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_KERNELS_INTERFACE_HPP
