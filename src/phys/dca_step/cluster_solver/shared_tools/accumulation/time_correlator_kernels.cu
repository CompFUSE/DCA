// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Implementation of the G0 computation for time measurements.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/kernels_interface.hpp"

#include "dca/util/cuda_blocks.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/solver_helper.cuh"
#include "dca/util/type_help.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace details {
// dca::phys::solver::details::

using namespace dca::linalg;
  using dca::util::SignType;
using dca::linalg::util::GpuStream;

template <typename Scalar, typename Real>
__global__ void computeG0Kernel(linalg::MatrixView<Scalar, linalg::GPU>& mat,
                                const DeviceInterpolationData<Scalar, Real> g0, const Real* t_l,
                                const int* b_l, const int* r_l, const Real* t_r, const int* b_r,
                                const int* r_r) {
  const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
  const unsigned j = blockIdx.y * blockDim.y + threadIdx.y;
  if (i >= mat.nrRows() || j >= mat.nrCols())
    return;

  const auto index = solver_helper.index(b_l[i], b_r[j], r_l[i], r_r[j]);
  const Real tau = t_l[i] - t_r[j];
  auto mat_elem_ptr = castGPUType(&mat(i, j));
  *mat_elem_ptr = castGPUType(g0(tau, index));
}

template <typename Scalar, typename Real, typename SignType>
void computeG0(linalg::MatrixView<Scalar, linalg::GPU>& g0_mat,
               const DeviceInterpolationData<Scalar, SignType> g0, const Real* t_l,
               const int* b_l, const int* r_l, const Real* t_r, const int* b_r, const int* r_r,
               const GpuStream& stream) {
  assert(SolverHelper::initialized());
  auto blocks = dca::util::get2DBlockSize(g0_mat.nrRows(), g0_mat.nrCols(), 32);

  computeG0Kernel<<<blocks[0], blocks[1], 0, stream>>>(g0_mat, g0, t_l, b_l, r_l, t_r, b_r, r_r);
}

// Instantation.
  template<> void computeG0<double, double, std::int8_t>(linalg::MatrixView<double, linalg::GPU>&,
					    const DeviceInterpolationData<double, std::int8_t>,
                                        const double*, const int*, const int*, const double*,
                                        const int*, const int*, const GpuStream&);
  template<> void computeG0<float, float, std::int8_t>(linalg::MatrixView<float, linalg::GPU>&,
					  const DeviceInterpolationData<float, std::int8_t>, const float*,
                                      const int*, const int*, const float*, const int*, const int*,
                                      const GpuStream&);
  template<> void computeG0<std::complex<double>, double, std::complex<double>>(
    linalg::MatrixView<std::complex<double>, linalg::GPU>&,
    const DeviceInterpolationData<std::complex<double>, std::complex<double>>, const double*, const int*,
    const int*, const double*, const int*, const int*, const GpuStream&);
  template<> void computeG0<std::complex<float>, float, std::complex<float>>(
    linalg::MatrixView<std::complex<float>, linalg::GPU>&,
    const DeviceInterpolationData<std::complex<float>, std::complex<float>>, const float*, const int*, const int*,
    const float*, const int*, const int*, const GpuStream&);


}  // namespace details
}  // namespace solver
}  // namespace phys
}  // namespace dca
