// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
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
// dca::phys::solver::ctint::

class DeviceInterpolationData {
public:
  DeviceInterpolationData(const DeviceInterpolationData& other) = default;

  __DEVICE__ double operator()(double tau, int lindex) const;

protected:
  DeviceInterpolationData() = default;

  constexpr static int coeff_size_ = 4;  // Same as G0Interpolation<linalg::CPU>.
  unsigned stride_;
  double beta_, n_div_beta_;
  double *values_, *g0_minus_;
};

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_DEVICE_INTERPOLATION_DATA_HPP
