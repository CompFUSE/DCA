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

// This amount of namespace nesting is really uncalled for
// http://open-std.org/JTC1/SC22/WG21/docs/papers/2017/p0816r0.pdf
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

// Updates G4 in the range [start, end)
template <typename Real, FourPointType type>
float updateG4(std::complex<Real>* G4, const std::complex<Real>* G_up, const int ldgu,
               const std::complex<Real>* G_down, const int ldgd, const int sign, bool atomic,
               cudaStream_t stream, std::size_t start, std::size_t end);

}  // namespace details
}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_KERNELS_INTERFACE_HPP
