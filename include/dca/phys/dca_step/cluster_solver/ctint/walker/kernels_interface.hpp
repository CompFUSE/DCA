// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_KERNELS_INTERFACE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_KERNELS_INTERFACE_HPP
#ifdef DCA_HAVE_CUDA

#include <cuda.h>

#include "dca/linalg/matrix_view.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

using MatrixView = linalg::MatrixView<double, linalg::GPU>;

void smallInverse(const MatrixView& m_in, MatrixView& m_out, cudaStream_t stream);

void smallInverse(MatrixView& in_out, cudaStream_t stream);

double separateIndexDeterminant(MatrixView m, const ushort* indices, int n_indices, double* det,
                         cudaStream_t stream);

}  // details
}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_KERNELS_INTERFACE_HPP
