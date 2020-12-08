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

#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/util/cuda_blocks.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/solver_helper.cuh"

namespace dca {
namespace phys {
namespace solver {
namespace details {
// dca::phys::solver::details::

using namespace dca::linalg;

template <typename Scalar, typename Real>
__global__ void computeG0Kernel(linalg::MatrixView<Scalar, linalg::GPU> mat,
                                const DeviceInterpolationData<Scalar> g0, const Real* t_l,
                                const int* b_l, const int* r_l, const Real* t_r, const int* b_r,
                                const int* r_r) {
  const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
  const unsigned j = blockIdx.y * blockDim.y + threadIdx.y;
  if (i >= mat.nrRows() || j >= mat.nrCols())
    return;

  const auto index = solver_helper.index(b_l[i], b_r[j], r_l[i], r_r[j]);
  const Real tau = t_l[i] - t_r[j];

  castCuda(mat(i, j)) = castCuda(g0(tau, index));
}

template <class T>
struct CudaScalarImpl {
  using type = T;
};
template <>
struct CudaScalarImpl<std::complex<double>> {
  using type = cuDoubleComplex;
};
template <>
struct CudaScalarImpl<std::complex<float>> {
  using type = cuComplex;
};

template <typename Scalar, typename Real>
void computeG0(linalg::MatrixView<Scalar, linalg::GPU>& g0_mat,
               const DeviceInterpolationData<Scalar> g0, const Real* t_l, const int* b_l,
               const int* r_l, const Real* t_r, const int* b_r, const int* r_r, cudaStream_t stream) {
  assert(SolverHelper::initialized());
  auto blocks = dca::util::get2DBlockSize(g0_mat.nrRows(), g0_mat.nrCols(), 32);
  using dca::linalg::castCuda;

  computeG0Kernel<<<blocks[0], blocks[1], 0, stream>>>(g0_mat, g0, t_l, b_l, r_l, t_r, b_r, r_r);
}

// Instantation.
template void computeG0(linalg::MatrixView<double, linalg::GPU>&,
                        const DeviceInterpolationData<double>, const double*, const int*,
                        const int*, const double*, const int*, const int*, cudaStream_t);
template void computeG0(linalg::MatrixView<float, linalg::GPU>&,
                        const DeviceInterpolationData<float>, const float*, const int*, const int*,
                        const float*, const int*, const int*, cudaStream_t);
template void computeG0(linalg::MatrixView<std::complex<double>, linalg::GPU>&,
                        const DeviceInterpolationData<std::complex<double>>, const double*,
                        const int*, const int*, const double*, const int*, const int*, cudaStream_t);
template void computeG0(linalg::MatrixView<std::complex<float>, linalg::GPU>&,
                        const DeviceInterpolationData<std::complex<float>>, const float*,
                        const int*, const int*, const float*, const int*, const int*, cudaStream_t);

}  // namespace details
}  // namespace solver
}  // namespace phys
}  // namespace dca
