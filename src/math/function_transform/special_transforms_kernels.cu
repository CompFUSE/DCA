// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Implements the GPU kernels used by the SpaceTransform2DGpu class.

#include "dca/math/function_transform/special_transforms/kernels_interface.hpp"

#include <array>

#include "dca/util/integer_division.hpp"
#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"

namespace dca {
namespace math {
namespace transform {
namespace details {
// dca::math::transform::details::

using linalg::util::castCudaComplex;
using linalg::util::CudaComplex;

std::array<dim3, 2> getBlockSize(const int i, const int j) {
  const int n_threads_i = std::min(32, i);
  const int n_threads_j = std::min(32, j);
  const int n_blocks_i = util::ceilDiv(i, n_threads_i);
  const int n_blocks_j = util::ceilDiv(j, n_threads_j);

  return std::array<dim3, 2>{dim3(n_blocks_i, n_blocks_j), dim3(n_threads_i, n_threads_j)};
}

template <typename Real>
__global__ void phaseFactorsAndRearrangeKernel(const CudaComplex<Real>* in, const int ldi,
                                               CudaComplex<Real>* out, const int ldo, const int nb,
                                               const int nk, const int nw,
                                               const CudaComplex<Real>* phase_factors) {
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;

  const int cols = nb * nk * nw;
  const int rows = cols / 2;
  if (id_i >= rows || id_j >= cols)
    return;

  const int no = nb * nk;
  auto get_indices = [nk, no](int id, int& b, int& k, int& w) {
    w = id / no;
    id -= w * no;
    b = id / nk;
    k = id - b * nk;
  };
  int w1, w2, b1, b2, k1, k2;

  get_indices(id_i, b1, k1, w1);
  get_indices(id_j, b2, k2, w2);

  const int out_i = b1 + nb * k1 + no * w1;
  const int out_j = b2 + nb * k2 + no * w2;

  CudaComplex<Real> value = in[id_i + ldi * id_j];

  using namespace dca::linalg;
  if (phase_factors)
    value *= phase_factors[b1 + nb * k1] * linalg::conj(phase_factors[b2 + nb * k2]);

  out[out_i + ldo * out_j] = value;
}

template <typename Real>
void phaseFactorsAndRearrange(const std::complex<Real>* in, const int ldi, std::complex<Real>* out,
                              const int ldo, const int nb, const int nk, const int nw,
                              const std::complex<Real>* phase_factors, const cudaStream_t stream) {
  const int size = nk * nb * nw;
  auto const blocks = getBlockSize(size / 2, size);

  phaseFactorsAndRearrangeKernel<Real>
      <<<blocks[0], blocks[1], 0, stream>>>(castCudaComplex(in), ldi, castCudaComplex(out), ldo, nb,
                                            nk, nw, castCudaComplex(phase_factors));
}

// Explicit instantiation.
template void phaseFactorsAndRearrange<double>(const std::complex<double>* in, const int ldi,
                                               std::complex<double>* out, const int ldo,
                                               const int nb, const int nk, const int nw,
                                               const std::complex<double>* phase_factors,
                                               const cudaStream_t stream);
template void phaseFactorsAndRearrange<float>(const std::complex<float>* in, const int ldi,
                                              std::complex<float>* out, const int ldo, const int nb,
                                              const int nk, const int nw,
                                              const std::complex<float>* phase_factors,
                                              cudaStream_t const stream);

}  // namespace details
}  // namespace transform
}  // namespace math
}  // namespace dca
