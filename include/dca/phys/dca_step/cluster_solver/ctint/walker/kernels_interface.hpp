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

void setRightSectorToId(double* m, int ldm, int n0, int n_max, cudaStream_t stream);

using MatrixView = linalg::MatrixView<double, linalg::GPU>;

void computeGLeft(MatrixView& G, const MatrixView& M, const double* f, int n_init,
                  cudaStream_t stream);

void multiplyByFFactor(MatrixView& M, const double* f_vals, bool inverse_factor, bool row_factor,
                       cudaStream_t stream);

void divideByGammaFactor(MatrixView m, const std::pair<int, double>* gamma_indices, int n_indices,
                         cudaStream_t stream);

}  // details
}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_KERNELS_INTERFACE_HPP
