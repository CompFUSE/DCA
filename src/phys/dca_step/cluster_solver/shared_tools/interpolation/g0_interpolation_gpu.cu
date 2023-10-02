// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the device methods of G0Interpolation<GPU>.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/kernels_interface.hpp"

#include "dca/platform/dca_gpu.h"

#include "dca/util/cuda_blocks.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/device_interpolation_data.hpp"
#include "dca/linalg/util/gpu_type_mapping.hpp"
namespace dca {
namespace phys {
namespace solver {
namespace details {

// dca::phys::solver::ctint::details::
using dca::util::CUDATypeMap;
using dca::util::castGPUType;
using dca::util::ComplexAlias;
using dca::util::RealAlias;
using dca::util::GPUTypeConversion;
using dca::util::IsComplex_t;
using dca::util::IsComplex;
using dca::util::IsReal;
using namespace dca::linalg;
using dca::util::SignType;
using dca::util::CudaComplex;
using dca::util::CudaScalar;

template <typename Scalar, typename Real, typename SignType>
__global__ void interpolateSlowKernel(Real tau, const int lindex,
                                      DeviceInterpolationData<Scalar, SignType> g0, CudaScalar<Scalar>* result) {
  *result = g0(tau, lindex);
}

template <typename Scalar, typename Real, typename SignType>
Scalar interpolateSlow(Real tau, int lindex,
                                          const DeviceInterpolationData<Scalar, SignType>& g0) {
  CudaScalar<Scalar>* d_result;
  Scalar result;
  cudaMalloc((void**)&d_result, sizeof(Scalar));

  interpolateSlowKernel<<<1, 1>>>(tau, lindex, g0, d_result);

  assert(cudaSuccess == cudaPeekAtLastError());
  cudaMemcpy(&result, d_result, sizeof(Scalar), cudaMemcpyDeviceToHost);
  cudaFree(d_result);
  return result;
}

template float interpolateSlow(float, int, const DeviceInterpolationData<float, std::int8_t>&);
template double interpolateSlow(double, int, const DeviceInterpolationData<double, std::int8_t>&);
  template std::complex<float> interpolateSlow(
    float, int, const DeviceInterpolationData<std::complex<float>, std::complex<float>>&);
  template std::complex<double> interpolateSlow(
    double, int, const DeviceInterpolationData<std::complex<double>, std::complex<double>>&);

}  // namespace details
}  // namespace solver
}  // namespace phys
}  // namespace dca
