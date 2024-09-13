// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//          Peter W. Doak (doakpw@ornl.gov)
//
// This file implements the device methods of DMatrixBuilder.

#include "dca/config/haves_defines.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/kernels_interface.hpp"
#include "dca/platform/dca_gpu.h"
#include "dca/platform/gpu_definitions.h"
#include "dca/linalg/util/gpu_type_mapping.hpp"
#include "dca/util/type_help.hpp"
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#include "dca/linalg/util/gpu_stream.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/solver_helper.cuh"
#include "dca/util/cuda_blocks.hpp"
#include <complex>

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

using dca::util::castGPUType;
using dca::linalg::assign;
using dca::util::SignType;
using dca::linalg::util::GpuStream;
  
template <typename Scalar, typename Real>
__global__ void buildG0MatrixKernel(linalg::MatrixView<Scalar, linalg::GPU> G0, const int n_init,
                                    const bool right_section, DeviceConfiguration config,
                                    DeviceInterpolationData<Scalar, Real> g0_interp) {
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;
  if (id_i >= G0.nrRows() || id_j >= G0.nrCols())
    return;
  int i(id_i);
  int j(id_j);

  if (right_section)
    j += n_init;
  else
    i += n_init;

  const int b_i = config.getLeftB(i);
  const auto tau_i = config.getTau(i);

  const int b_j = config.getRightB(j);
  const auto tau_j = config.getTau(j);

  const int label = dca::phys::solver::details::solver_helper.index(b_i, b_j, config.getLeftR(i),
                                                                    config.getRightR(j));

  auto g0_elem_ptr = castGPUType(&G0(id_i, id_j));
  assign(*g0_elem_ptr, g0_interp(tau_i - tau_j, label));
}

template <typename Scalar, typename SignType>
void buildG0Matrix(linalg::MatrixView<Scalar, linalg::GPU> g0, const int n_init,
                   const bool right_section, DeviceConfiguration config,
                   DeviceInterpolationData<Scalar, SignType> g0_interp, const GpuStream& stream) {
  // assert(CtintHelper::is_initialized());
  const auto blocks = dca::util::getBlockSize(g0.nrRows(), g0.nrCols());

  buildG0MatrixKernel<<<blocks[0], blocks[1], 0, stream.streamActually()>>>(g0, n_init, right_section, config,
                                                           g0_interp);
  checkErrorsCudaDebug();
}

template void buildG0Matrix(linalg::MatrixView<float, linalg::GPU>, const int, const bool,
                            DeviceConfiguration, DeviceInterpolationData<float, signed char>, const GpuStream&);
template void buildG0Matrix(linalg::MatrixView<double, linalg::GPU>, const int, const bool,
                            DeviceConfiguration, DeviceInterpolationData<double, std::int8_t>, const GpuStream&);
template void buildG0Matrix(linalg::MatrixView<std::complex<float>, linalg::GPU>, const int,
                            const bool, DeviceConfiguration,
                            DeviceInterpolationData<std::complex<float>, std::complex<float>>, const GpuStream&);
template void buildG0Matrix(linalg::MatrixView<std::complex<double>, linalg::GPU>, const int,
                            const bool, DeviceConfiguration,
                            DeviceInterpolationData<std::complex<double>, std::complex<double>>, const GpuStream&);

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
