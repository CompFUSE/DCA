// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a trivially copyable representation of the data needed by
// g0_interpolation_gpu.hpp.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_DEVICE_INTERPOLATION_DATA_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_DEVICE_INTERPOLATION_DATA_HPP
#ifdef DCA_HAVE_CUDA

#include "dca/util/cuda_definitions.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

class DeviceInterpolationData {
public:
  DeviceInterpolationData(const DeviceInterpolationData& other) = default;

  __DEVICE__ double operator()(double tau, int lindex) const;


protected:
  DeviceInterpolationData() = default;

  constexpr static int coeff_size_ = 4; // Same as G0Interpolation<linalg::CPU>.
  double beta_, n_div_beta_;
  double *values_, *g0_minus_;
};

}  // dca
}  // phys
}  // solver
}  // ctint
}  // details

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_DEVICE_INTERPOLATION_DATA_HPP
