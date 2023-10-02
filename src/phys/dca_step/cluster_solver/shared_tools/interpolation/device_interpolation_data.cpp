// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//          Peter Doak (doakpw@ornl.gov)
//
// This file provides explicit instantiation for trivially copyable representation of the data needed by
// g0_interpolation_gpu.hpp.


#include "dca/linalg/util/complex_operators_cuda.cu.hpp"


namespace dca {
namespace phys {
namespace solver {

template<typename Scalar, typename Sign>
__DEVICE__ auto DeviceInterpolationData::operator()(Real tau, int lindex) const {
  assert(tau >= -beta_ && tau <= beta_);

  if (tau == 0)  // returns G0(tau = 0+)
    return *castGPUType(g0_minus_ + lindex);

  Sign factor{};
  
  if (tau < 0) {
    tau += beta_;
    factor += -1.0;
  }
  else
    factor += 1.0;

  // Scale tau in [0, n_time_slices). Assume even spacing in time.
  const Real scaled_tau = tau * n_div_beta_;
  const int tau_index(scaled_tau);
  const Real delta_tau = scaled_tau - tau_index;

  // Get the pointer to the first akima coeff.
  const CUDATypeMap<Scalar>* coeff_ptr =
      castGPUType(values_) + (tau_index * coeff_size_ + lindex * stride_);

  // Return akima interpolation.
  return factor *
    static_cast<decltype<factor>>(coeff_ptr[0] +
          delta_tau * (coeff_ptr[1] + delta_tau * (coeff_ptr[2] + delta_tau * coeff_ptr[3])));
}

 template class DeviceInterpolationData<float, int>;
 template class DeviceInterpolationData<double, int>;
 template class DeviceInterpolationData<std::complex<float>, std : complex<float>>;
 template class DeviceInterpolationData<std::complex<double>, std : complex<double>>;

 template std::complex<double> DeviceInterpolationData<std::complex<double>, std::complex<double>>::operator()(double, int);
 template std::complex<float> DeviceInterpolationData<std::complex<float>, std::complex<float>>::operator()(float, int);
 template float DeviceInterpolationData<float, float>::operator()(float, int);
 template double DeviceInterpolationData<double, double>::operator()(double, int);

}  // namespace solver
}  // namespace phys
}  // namespace dca
