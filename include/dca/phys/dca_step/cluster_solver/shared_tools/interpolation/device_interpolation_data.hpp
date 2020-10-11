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

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_DEVICE_INTERPOLATION_DATA_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_DEVICE_INTERPOLATION_DATA_HPP
#ifdef DCA_HAVE_CUDA

#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#include "dca/util/cuda_definitions.hpp"
#include "dca/util/type_utils.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

using linalg::CudaConvert;

template <typename Scalar>
class DeviceInterpolationData {
public:
  using Real = dca::util::Real<Scalar>;

  DeviceInterpolationData(const DeviceInterpolationData& other) = default;

  __DEVICE__ auto operator()(Real tau, int lindex) const {
    using namespace dca::linalg;

    assert(tau >= -beta_ && tau <= beta_);

    if (tau == 0)  // returns G0(tau = 0+)
      return castCuda(g0_minus_[lindex]);

    Real factor = 1;
    if (tau < 0) {
      tau += beta_;
      factor = -1;
    }

    // Scale tau in [0, n_time_slices). Assume even spacing in time.
    const Real scaled_tau = tau * n_div_beta_;
    const int tau_index(scaled_tau);
    const Real delta_tau = scaled_tau - tau_index;

    // Get the pointer to the first akima coeff.
    const CudaConvert<Scalar>* coeff_ptr =
        castCuda(values_) + (tau_index * coeff_size_ + lindex * stride_);

    // Return akima interpolation.
    return (coeff_ptr[0] +
            delta_tau * (coeff_ptr[1] + delta_tau * (coeff_ptr[2] + delta_tau * coeff_ptr[3]))) *
           factor;
  }

protected:
  DeviceInterpolationData() = default;

  constexpr static int coeff_size_ = 4;  // Same as G0Interpolation<linalg::CPU>.
  unsigned stride_;
  Real beta_, n_div_beta_;
  Scalar* values_;
  Scalar* g0_minus_;
};

}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_DEVICE_INTERPOLATION_DATA_HPP
