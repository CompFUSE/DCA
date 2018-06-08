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

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_INTERFACE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_INTERFACE_HPP
#ifdef DCA_HAVE_CUDA

#include <cuda.h>

#include "dca/linalg/matrix_view.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/device_interpolation_data.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

using Configuration = DeviceConfiguration;
using Interpolation = details::DeviceInterpolationData;
using MatrixView = linalg::MatrixView<double, linalg::GPU>;

void buildG0Matrix(MatrixView G0, int n_init, bool right_section, Configuration config,
                   Interpolation g0_interp, cudaStream_t stream);

double deviceInterpolationTest(DeviceInterpolationData data, double tau, int linindex);

}  // details
}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_INTERFACE_HPP
