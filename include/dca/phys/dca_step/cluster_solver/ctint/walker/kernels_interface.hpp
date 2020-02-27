// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_KERNELS_INTERFACE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_KERNELS_INTERFACE_HPP

#include <cuda.h>

#include "dca/linalg/matrix_view.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {

template <typename Real>
void setRightSectorToId(Real* m, int ldm, int n0, int n_max, cudaStream_t stream);

template <typename Real>
using MatrixView = linalg::MatrixView<Real, linalg::GPU>;

template <typename Real>
void computeGLeft(MatrixView<Real>& G, const MatrixView<Real>& M, const Real* f, int n_init,
                  cudaStream_t stream);

template <typename Real>
void multiplyByFColFactor(MatrixView<Real>& M, const Real* f_vals, cudaStream_t stream);

template <typename Real>
void multiplyByInverseFFactor(const MatrixView<Real>& m_in, MatrixView<Real>& m_out,
                              const Real* f_vals, cudaStream_t stream);

template <typename Real>
void divideByGammaFactor(MatrixView<Real> m, const std::pair<int, Real>* gamma_indices,
                         int n_indices, cudaStream_t stream);

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_KERNELS_INTERFACE_HPP
