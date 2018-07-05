// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Implements the GPU kernels used by 'cached_ndft_gpu.hpp'.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/kernels_interface.hpp"

#include <array>
#include <cassert>
#include <cuda.h>
#include <cuda_runtime.h>

#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/util/cuda_blocks.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {

using linalg::util::CudaComplex;
using linalg::util::castCudaComplex;

template <typename Real, typename InpScalar>
__global__ void sortMKernel(const int size, const InpScalar* M, const int ldm,
                            CudaComplex<Real>* sorted_M, int lds, const Triple<Real>* config1,
                            const Triple<Real>* config2) {
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;
  if (id_i >= size || id_j >= size)
    return;

  const int inp_i = config1[id_i].idx;
  const int inp_j = config2[id_j].idx;

  sorted_M[id_i + lds * id_j].x = M[inp_i + ldm * inp_j];
  sorted_M[id_i + lds * id_j].y = 0;
}

template <typename Real, typename InpScalar>
void sortM(const int size, const InpScalar* M, const int ldm, std::complex<Real>* sorted_M,
           const int lds, const Triple<Real>* config1, const Triple<Real>* config2,
           const cudaStream_t stream) {
  if (!size)
    return;

  auto const blocks = dca::util::getBlockSize(size, size);

  sortMKernel<<<blocks[0], blocks[1], 0, stream>>>(size, M, ldm, castCudaComplex(sorted_M), lds,
                                                   config1, config2);
}

template <typename Real>
__global__ void computeTKernel(const int n, const int m, CudaComplex<Real>* T, int ldt,
                               const Triple<Real>* config, const Real* w, const bool transposed) {
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;
  if (id_i >= n || id_j >= m)
    return;

  if (!transposed) {
    T[id_i + ldt * id_j].x = cos(w[id_i] * config[id_j].tau);
    T[id_i + ldt * id_j].y = sin(w[id_i] * config[id_j].tau);
  }
  else {
    T[id_i + ldt * id_j].x = cos(w[id_j] * config[id_i].tau);
    T[id_i + ldt * id_j].y = -sin(w[id_j] * config[id_i].tau);
  }
}

template <typename Real>
void computeT(const int n, const int m, std::complex<Real>* T, int ldt, const Triple<Real>* config,
              const Real* w, const bool transposed, const cudaStream_t stream) {
  auto const blocks = dca::util::getBlockSize(n, m);

  computeTKernel<<<blocks[0], blocks[1], 0, stream>>>(n, m, castCudaComplex(T), ldt, config, w,
                                                      transposed);
}

template <typename Real>
__global__ void rearrangeOutputKernel(const int nw, const int no, const int nb,
                                      const CudaComplex<Real>* in, const int ldi,
                                      CudaComplex<Real>* out, const int ldo) {
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;
  const int n_rows = nw / 2 * no;
  const int n_cols = nw * no;
  if (id_i >= n_rows || id_j >= n_cols)
    return;

  auto get_indices = [nb](int id, const int nw, int& b, int& r, int& w) {
    r = id / (nw * nb);
    id -= r * nw * nb;
    b = id / nw;
    w = id - b * nw;
  };
  int w1, w2, b1, b2, r1, r2;

  get_indices(id_i, nw / 2, b1, r1, w1);
  get_indices(id_j, nw, b2, r2, w2);

  const int nr = no / nb;
  const int out_i = r1 + nr * b1 + no * w1;
  const int out_j = r2 + nr * b2 + no * w2;

  out[out_i + ldo * out_j] = in[id_i + ldi * id_j];
}

template <typename Real>
void rearrangeOutput(const int nw, const int no, const int nb, const std::complex<Real>* in,
                     const int ldi, std::complex<Real>* out, const int ldo,
                     const cudaStream_t stream) {
  const int n_rows = nw / 2 * no;
  const int n_cols = nw * no;
  auto const blocks = dca::util::getBlockSize(n_rows, n_cols);

  rearrangeOutputKernel<Real><<<blocks[0], blocks[1], 0, stream>>>(nw, no, nb, castCudaComplex(in),
                                                                   ldi, castCudaComplex(out), ldo);
}

// Explicit instantiation.
template void sortM<double, double>(int, const double*, int, std::complex<double>*, const int,
                                    const Triple<double>*, const Triple<double>*, const cudaStream_t);
template void sortM<float, double>(int, const double*, int, std::complex<float>*, const int,
                                   const Triple<float>*, const Triple<float>*, const cudaStream_t);
template void computeT<double>(int, int, std::complex<double>*, int, const Triple<double>*,
                               const double*, bool, const cudaStream_t);
template void computeT<float>(int, int, std::complex<float>*, int, const Triple<float>*,
                              const float*, bool, const cudaStream_t);
template void rearrangeOutput<double>(const int nw, const int no, const int nb,
                                      const std::complex<double>* in, const int ldi,
                                      std::complex<double>* out, const int ldo,
                                      const cudaStream_t stream);
template void rearrangeOutput<float>(const int nw, const int no, const int nb,
                                     const std::complex<float>* in, const int ldi,
                                     std::complex<float>* out, const int ldo,
                                     const cudaStream_t stream);

}  // details
}  // accumulator
}  // solver
}  // phys
}  // dca
