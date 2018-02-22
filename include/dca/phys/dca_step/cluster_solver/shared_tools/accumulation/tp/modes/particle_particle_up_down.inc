// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Specialized kernels for the PARTICL_PARTICLE_UP_DOWN_MODE
// G4 definition: G4 = 1/2 \sum_s <c^+(q-k1, s) c^+(k1, -s) c(k2, -s) c(q-k2, s)>

template <typename Real>
__global__ void updateG4PPUpDownSinglebandKernel(CudaComplex<Real>* __restrict__ G4,
                                                 const CudaComplex<Real>* __restrict__ G_up,
                                                 const int ldgu,
                                                 const CudaComplex<Real>* __restrict__ G_down,
                                                 const int ldgd, const int nk, const int nw,
                                                 const int sign, const G4Helper helper) {
  const int size = nk * nw;
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;
  if (id_i >= size || id_j >= size)
    return;

  auto get_indices = [=](const int id, int& k, int& w) {
    w = id / nk;
    k = id - nk * w;
  };
  int w1, w2, k1, k2;
  get_indices(id_i, k1, w1);
  get_indices(id_j, k2, w2);

  int w1_a(w1), w2_a(w2);
  const bool conj_a = helper.extendWIndices(w1_a, w2_a);
  const int i_a = k1 + nk * w1_a;
  const int j_a = k2 + nk * w2_a;

  int w1_b(helper.wMinus(w1));
  int w2_b(helper.wMinus(w2));
  const bool conj_b = helper.extendWIndices(w1_b, w2_b);
  const int i_b = helper.qMinus(k1) + nk * w1_b;
  const int j_b = helper.qMinus(k2) + nk * w2_b;

  auto cond_conj = [](const CudaComplex<Real> a, const bool cond) { return cond ? conj(a) : a; };
  const Real factor = 0.5 * sign;

  const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
  const CudaComplex<Real> Gb_1 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);

  const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);
  const CudaComplex<Real> Gb_2 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

  G4[helper.g4Index(k1, k2, w1, w2)] += (Ga_1 * Gb_1 + Ga_2 * Gb_2) * factor;
}

template <typename Real>
__global__ void updateG4PPUpDownMultibandKernel(CudaComplex<Real>* __restrict__ G4,
                                                const CudaComplex<Real>* __restrict__ G_up,
                                                const int ldgu,
                                                const CudaComplex<Real>* __restrict__ G_down,
                                                const int ldgd, const int nb, const int nk,
                                                const int nw, const int sign, const G4Helper helper) {
  const int size = nk * nw * nb * nb;
  const int id_i = blockIdx.x * blockDim.x + threadIdx.x;
  const int id_j = blockIdx.y * blockDim.y + threadIdx.y;
  if (id_i >= size || id_j >= size)
    return;

  const int step2 = nb * nb;
  const int step1 = step2 * nk;
  auto get_indices = [=](int id, int& b1, int& b2, int& k, int& w) {
    w = id / step1;
    id -= w * step1;
    k = id / step2;
    id -= k * step2;
    b2 = id / nb;
    b1 = id - nb * b2;
  };
  int w1, w2, k1, k2, b1, b2, b3, b4;
  get_indices(id_i, b1, b2, k1, w1);
  get_indices(id_j, b3, b4, k2, w2);

  int w1_a(w1), w2_a(w2);
  const int no = nk * nb;
  const bool conj_a = helper.extendWIndices(w1_a, w2_a);
  const int i_a = b1 + nb * k1 + no * w1_a;
  const int j_a = b2 + nb * k2 + no * w2_a;

  int w1_b(helper.wMinus(w1));
  int w2_b(helper.wMinus(w2));
  const bool conj_b = helper.extendWIndices(w1_b, w2_b);
  const int i_b = b3 + nb * helper.qMinus(k1) + no * w1_b;
  const int j_b = b4 + nb * helper.qMinus(k2) + no * w2_b;

  auto cond_conj = [](const CudaComplex<Real> a, const bool cond) { return cond ? conj(a) : a; };
  const Real factor = 0.5 * sign;

  const CudaComplex<Real> Ga_1 = cond_conj(G_up[i_a + ldgu * j_a], conj_a);
  const CudaComplex<Real> Gb_1 = cond_conj(G_down[i_b + ldgd * j_b], conj_b);

  const CudaComplex<Real> Ga_2 = cond_conj(G_down[i_a + ldgd * j_a], conj_a);
  const CudaComplex<Real> Gb_2 = cond_conj(G_up[i_b + ldgu * j_b], conj_b);

  G4[helper.g4Index(b1, b3, b2, b4, k1, k2, w1, w2)] += (Ga_1 * Gb_1 + Ga_2 * Gb_2) * factor;
}

template <typename Real>
void updateG4PPUpDown(std::complex<Real>* G4, const std::complex<Real>* G_up, const int ldgu,
                      const std::complex<Real>* G_down, const int ldgd, const int nb, const int nk,
                      const int nw_pos, const int sign, cudaStream_t stream) {
  const int nw = 2 * nw_pos;
  const int size = nw * nk * nb * nb;
  const static auto blocks = getBlockSize(size, size, 16);

  if (nb == 1)
    updateG4PPUpDownSinglebandKernel<Real><<<blocks[0], blocks[1], 0, stream>>>(
        castCudaComplex(G4), castCudaComplex(G_up), ldgu, castCudaComplex(G_down), ldgd, nk, nw,
        sign, global::helper);
  else
    updateG4PPUpDownMultibandKernel<Real><<<blocks[0], blocks[1], 0, stream>>>(
        castCudaComplex(G4), castCudaComplex(G_up), ldgu, castCudaComplex(G_down), ldgd, nb, nk, nw,
        sign, global::helper);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

// Explicit instantiation.
template void updateG4PPUpDown<float>(std::complex<float>* G4, const std::complex<float>* G_up,
                                      const int ldgu, const std::complex<float>* G_down,
                                      const int ldgd, const int nb, const int nk, const int nw_pos,
                                      const int sign, cudaStream_t stream);
template void updateG4PPUpDown<double>(std::complex<double>* G4, const std::complex<double>* G_up,
                                       const int ldgu, const std::complex<double>* G_down,
                                       const int ldgd, const int nb, const int nk, const int nw_pos,
                                       const int sign, cudaStream_t stream);
