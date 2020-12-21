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

template <typename Real>
class DeviceInterpolationData {
public:
  DeviceInterpolationData(const DeviceInterpolationData& other) = default;

  __DEVICE__ Real operator()(Real tau, int lindex) const {
    assert(tau >= -beta_ && tau <= beta_);

    if (tau == 0)  // returns G0(tau = 0+)
      return g0_minus_[lindex];

    short int factor = 1;
    if (tau < 0) {
      tau += beta_;
      factor = -1;
    }

    // Scale tau in [0, n_time_slices). Assume even spacing in time.
    const Real scaled_tau = tau * n_div_beta_;
    const int tau_index(scaled_tau);
    const Real delta_tau = scaled_tau - tau_index;

    // Get the pointer to the first akima coeff.
    const Real* coeff_ptr = &values_[tau_index * coeff_size_ + lindex * stride_];

    // Return akima interpolation.
    return factor *
           (coeff_ptr[0] +
            delta_tau * (coeff_ptr[1] + delta_tau * (coeff_ptr[2] + delta_tau * coeff_ptr[3])));
  }

protected:
  DeviceInterpolationData() = default;

  constexpr static int coeff_size_ = 4;  // Same as G0Interpolation<linalg::CPU>.
  unsigned stride_;
  Real beta_, n_div_beta_;
  Real *values_, *g0_minus_;
};

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_DEVICE_INTERPOLATION_DATA_HPP
