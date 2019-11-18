// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements g0_interpolation_kernels.hpp.

#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_kernels.hpp"

#include <cassert>

#include "cuda_runtime.h"

#include "dca/linalg/util/error_cuda.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/util/integer_division.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
namespace g0kernels {
// dca::phys::solver::ctaux::g0kernels::

__global__ void interpolation_kernel(int Nb, int Nr, int Nt, double beta, int Nv, int* b, int* r,
                                     double* t, double* G0, std::pair<int, int> G0_cs,
                                     std::pair<int, int> G0_gs, double* r0_min_r1,
                                     std::pair<int, int> r0_min_r1_cs,
                                     std::pair<int, int> r0_min_r1_gs, double* G0_r_t,
                                     std::pair<int, int> G0_r_t_cs, std::pair<int, int> G0_r_t_gs,
                                     double* grad_G0_r_t, std::pair<int, int> grad_G0_r_t_cs,
                                     std::pair<int, int> grad_G0_r_t_gs) {
  int I = threadIdx.x + blockDim.x * blockIdx.x;
  int J = threadIdx.y + blockDim.y * blockIdx.y;

  if (I > -1 && I < Nv && J > -1 && J < Nv) {
    int delta_r = r0_min_r1[r[J] + r[I] * r0_min_r1_gs.first];
    double tau = t[I] - t[J];

    double scaled_tau = (tau + beta) * double(Nt) / (2. * beta);

    int t_ind = scaled_tau;
    double delta_tau = scaled_tau - t_ind;

    int linind = b[I] + Nb * (b[J] + Nb * delta_r);

    double f_0 = G0_r_t[t_ind + linind * G0_r_t_gs.first];
    double grad = grad_G0_r_t[t_ind + linind * grad_G0_r_t_gs.first];

    G0[I + G0_gs.first * J] = -(f_0 + grad * delta_tau);

    if (I == J) {
      int t0_ind = Nt / 2;
      int r0_ind = r0_min_r1[0];

      G0[I + G0_gs.first * I] =
          -G0_r_t[t0_ind + G0_r_t_gs.first * (b[I] + Nb * (b[I] + Nb * r0_ind))];
    }
  }
}

void interpolate_G0_matrix_on_GPU(int Nb, int Nr, int Nt, double beta, int Nv, int* b, int* r,
                                  double* t, double* G0, std::pair<int, int> G0_cs,
                                  std::pair<int, int> G0_gs, double* r0_min_r1,
                                  std::pair<int, int> r0_min_r1_cs,
                                  std::pair<int, int> r0_min_r1_gs, double* G0_r_t,
                                  std::pair<int, int> G0_r_t_cs, std::pair<int, int> G0_r_t_gs,
                                  double* grad_G0_r_t, std::pair<int, int> grad_G0_r_t_cs,
                                  std::pair<int, int> grad_G0_r_t_gs) {
  checkErrorsCudaDebug();

  int Nr_t = 16;
  int Nr_b = dca::util::ceilDiv(Nv, Nr_t);

  dim3 threads(Nr_t, Nr_t);
  dim3 blocks(Nr_b, Nr_b);

  interpolation_kernel<<<blocks, threads>>>(Nb, Nr, Nt, beta, Nv, b, r, t, G0, G0_cs, G0_gs,
                                            r0_min_r1, r0_min_r1_cs, r0_min_r1_gs, G0_r_t, G0_r_t_cs,
                                            G0_r_t_gs, grad_G0_r_t, grad_G0_r_t_cs, grad_G0_r_t_gs);

  checkErrorsCudaDebug();
}

/***********************************
 ***
 ***        LINEAR-interpolation
 ***
 ************************************/

/*
  const static int BLOCK_SIZE        = 32;
  const static int NUMBER_OF_THREADS = 64;

  __global__ void interpolation_kernel(int Nb,
  int Nr,
  int Nt,
  double beta,
  int     Nc,
  int     Nv,
  int*    b,
  int*    r,
  double* t,
  double* G0         , std::pair<int,int> G0_cs         , std::pair<int,int> G0_gs,
  double* r0_min_r1  , std::pair<int,int> r0_min_r1_cs  , std::pair<int,int> r0_min_r1_gs,
  double* G0_r_t     , std::pair<int,int> G0_r_t_cs     , std::pair<int,int> G0_r_t_gs,
  double* grad_G0_r_t, std::pair<int,int> grad_G0_r_t_cs, std::pair<int,int> grad_G0_r_t_gs)
  {
  int I = threadIdx.x + blockDim.x*blockIdx.x;
  int J = threadIdx.y + blockDim.y*blockIdx.y;

  if( (I<Nv && J>=Nc && J<Nv) or
  (I>=Nc && I<Nv && J<Nv) )
  {
  int    delta_r = r0_min_r1[ r[J] + r[I]*r0_min_r1_gs.first ];
  double tau     = t[I]-t[J];

  double scaled_tau = (tau+beta)*double(Nt)/(2.*beta);

  int    t_ind     = scaled_tau;
  double delta_tau = scaled_tau-t_ind;

  int linind = b[I]+Nb*(b[J]+Nb*delta_r);

  double f_0  =      G0_r_t[t_ind+linind*     G0_r_t_gs.first];
  double grad = grad_G0_r_t[t_ind+linind*grad_G0_r_t_gs.first];

  G0[I+G0_gs.first*J] = -(f_0 + grad*delta_tau);

  if(I==J)
  {
  int t0_ind = Nt/2;
  int r0_ind = r0_min_r1[0];

  G0[I+G0_gs.first*I] = -G0_r_t[t0_ind+G0_r_t_gs.first*(b[I]+Nb*(b[I]+Nb*r0_ind))];
  }
  }
  }

  __global__ void interpolation_kernel_fat_column(int Nb,
  int Nr,
  int Nt,
  double beta,
  int     Nc,
  int     Nv,
  int*    b,
  int*    r,
  double* t,
  double* G0         , std::pair<int,int> G0_cs         , std::pair<int,int> G0_gs,
  double* r0_min_r1  , std::pair<int,int> r0_min_r1_cs  , std::pair<int,int> r0_min_r1_gs,
  double* G0_r_t     , std::pair<int,int> G0_r_t_cs     , std::pair<int,int> G0_r_t_gs,
  double* grad_G0_r_t, std::pair<int,int> grad_G0_r_t_cs, std::pair<int,int> grad_G0_r_t_gs)
  {
  int I     = threadIdx.x + blockDim.x*blockIdx.x;

  int J_min = Nc;
  int J_max = Nv;

  if(I>-1 && I<Nv)
  {
  for(int J=J_min; J<J_max; ++J)
  {
  int    delta_r = r0_min_r1[ r[J] + r[I]*r0_min_r1_gs.first ];
  double tau     = t[I]-t[J];

  double scaled_tau = (tau+beta)*double(Nt)/(2.*beta);

  int    t_ind     = scaled_tau;
  double delta_tau = scaled_tau-t_ind;

  int linind = b[I]+Nb*(b[J]+Nb*delta_r);

  double f_0  =      G0_r_t[t_ind+linind*     G0_r_t_gs.first];
  double grad = grad_G0_r_t[t_ind+linind*grad_G0_r_t_gs.first];

  G0[I+G0_gs.first*J] = -(f_0 + grad*delta_tau);


  //if(false and I==J)
  //G0[I+G0_gs.first*J] = -G0_r_t[t0_ind+G0_r_t_gs.first*(b[I]+Nb*(b[J]+Nb*r0_ind))];
  }
  }

  if(I>=J_min && I<J_max) // I==J
  {
  int t0_ind = Nt/2;
  int r0_ind = r0_min_r1[0];

  G0[I+G0_gs.first*I] = -G0_r_t[t0_ind+G0_r_t_gs.first*(b[I]+Nb*(b[I]+Nb*r0_ind))];
  }
  }

  __global__ void interpolation_kernel_fat_row(int Nb,
  int Nr,
  int Nt,
  double beta,
  int     Nc,
  int     Nv,
  int*    b,
  int*    r,
  double* t,
  double* G0         , std::pair<int,int> G0_cs         , std::pair<int,int> G0_gs,
  double* r0_min_r1  , std::pair<int,int> r0_min_r1_cs  , std::pair<int,int> r0_min_r1_gs,
  double* G0_r_t     , std::pair<int,int> G0_r_t_cs     , std::pair<int,int> G0_r_t_gs,
  double* grad_G0_r_t, std::pair<int,int> grad_G0_r_t_cs, std::pair<int,int> grad_G0_r_t_gs)
  {
  int I     = Nc+threadIdx.x;

  int J_min = BLOCK_SIZE*(blockIdx.x+0);
  int J_max = BLOCK_SIZE*(blockIdx.x+1);

  J_min = max(J_min, 0);
  J_max = min(J_max, Nc);

  if(I>=Nc && I<Nv)
  {
  for(int J=J_min; J<J_max; ++J)
  {
  int    delta_r = r0_min_r1[ r[J] + r[I]*r0_min_r1_gs.first ];
  double tau     = t[I]-t[J];

  double scaled_tau = (tau+beta)*double(Nt)/(2.*beta);

  int    t_ind     = scaled_tau;
  double delta_tau = scaled_tau-t_ind;

  int linind = b[I]+Nb*(b[J]+Nb*delta_r);

  double f_0  =      G0_r_t[t_ind+linind*     G0_r_t_gs.first];
  double grad = grad_G0_r_t[t_ind+linind*grad_G0_r_t_gs.first];

  G0[I+G0_gs.first*J] = -(f_0 + grad*delta_tau);
  }
  }
  }

  void interpolate_G0_matrix_on_GPU(int Nb,
  int Nr,
  int Nt,
  double beta,
  int     Nc,
  int     Nv,
  int*    b,
  int*    r,
  double* t,
  double* G0         , std::pair<int,int> G0_cs         , std::pair<int,int> G0_gs,
  double* r0_min_r1  , std::pair<int,int> r0_min_r1_cs  , std::pair<int,int> r0_min_r1_gs,
  double* G0_r_t     , std::pair<int,int> G0_r_t_cs     , std::pair<int,int> G0_r_t_gs,
  double* grad_G0_r_t, std::pair<int,int> grad_G0_r_t_cs, std::pair<int,int> grad_G0_r_t_gs)
  {
  //assert(cuda_check_for_errors("init 2 interpolation_kernel"));

  {
  int Nr_t = NUMBER_OF_THREADS;
  int Nr_b = dca::util::ceilDiv(Nv, Nr_t);

  dim3 threads(Nr_t);
  dim3 blocks (Nr_b);

  interpolation_kernel_fat_column<<<blocks, threads>>>(Nb, Nr, Nt, beta, Nc, Nv, b, r, t,
  G0         , G0_cs         , G0_gs         ,
  r0_min_r1  , r0_min_r1_cs  , r0_min_r1_gs  ,
  G0_r_t     , G0_r_t_cs     , G0_r_t_gs     ,
  grad_G0_r_t, grad_G0_r_t_cs, grad_G0_r_t_gs);

  //assert(cuda_check_for_errors("interpolation_kernel_fat_column"));
  }

  if(Nv-Nc>0)
  {
  int Nr_t = Nv-Nc;
  int Nr_b = dca::util::ceilDiv(Nc, BLOCK_SIZE);

  dim3 threads(Nr_t);
  dim3 blocks (Nr_b);

  interpolation_kernel_fat_row<<<blocks, threads>>>(Nb, Nr, Nt, beta, Nc, Nv, b, r, t,
  G0         , G0_cs         , G0_gs         ,
  r0_min_r1  , r0_min_r1_cs  , r0_min_r1_gs  ,
  G0_r_t     , G0_r_t_cs     , G0_r_t_gs     ,
  grad_G0_r_t, grad_G0_r_t_cs, grad_G0_r_t_gs);

  //assert(cuda_check_for_errors("interpolation_kernel_fat_row"));
  }


  //       {
  //      int Nr_t = 16;
  //      int Nr_b = dca::util::ceilDiv(Nv, Nr_t);

  //      dim3 threads(Nr_t, Nr_t);
  //      dim3 blocks (Nr_b, Nr_b);

  //      interpolation_kernel<<<blocks, threads>>>(Nb, Nr, Nt, beta, Nc, Nv, b, r, t,
  //                                                G0         , G0_cs         , G0_gs         ,
  //                                                r0_min_r1  , r0_min_r1_cs  , r0_min_r1_gs  ,
  //                                                G0_r_t     , G0_r_t_cs     , G0_r_t_gs     ,
  //                                                grad_G0_r_t, grad_G0_r_t_cs, grad_G0_r_t_gs);
  //       }

  }
*/

/***********************************
 ***
 ***        AKIMA-interpolation
 ***
 ************************************/

const static int BLOCK_SIZE_x = 32;
const static int BLOCK_SIZE_y = 16;

template<typename Real>
__global__ void akima_interpolation_fat_column(int Nb, int Nr, int Nt, Real beta, int Nc, int Nv,
                                               int* b, int* r, Real* t, Real* G0,
                                               std::pair<int, int> G0_cs, std::pair<int, int> G0_gs,
                                               Real* r0_min_r1, std::pair<int, int> r0_min_r1_cs,
                                               std::pair<int, int> r0_min_r1_gs, Real* alpha,
                                               std::pair<int, int> alpha_cs,
                                               std::pair<int, int> alpha_gs) {
  assert(blockDim.x == BLOCK_SIZE_x);

  int I = threadIdx.x + BLOCK_SIZE_x * blockIdx.x;

  int J_min = Nc + BLOCK_SIZE_y * (blockIdx.y + 0);
  int J_max = Nc + BLOCK_SIZE_y * (blockIdx.y + 1);

  J_min = max(Nc, J_min);
  J_max = min(Nv, J_max);

  if (I > -1 && I < Nv) {
    for (int J = J_min; J < J_max; ++J) {
      int delta_r = r0_min_r1[r[J] + r[I] * r0_min_r1_gs.first];
      Real tau = t[I] - t[J];

      Real scaled_tau = (tau + beta) * Real(Nt) / (2. * beta);

      int t_ind = scaled_tau;
      Real delta_tau = scaled_tau - t_ind;

      assert(delta_tau > -1.e-16 and delta_tau < 1 + -1.e-16);

      int col_ind = b[I] + Nb * (b[J] + Nb * delta_r);
      int row_ind = 4 * t_ind;
      int alpha_LD = alpha_gs.first;

      assert(row_ind > -1 and row_ind < alpha_cs.first);
      assert(col_ind > -1 and col_ind < alpha_cs.second);

      Real* a_ptr = &alpha[row_ind + col_ind * alpha_LD];

      assert(I > -1 and I < G0_cs.first);
      assert(J > -1 and J < G0_cs.second);

      G0[I + G0_gs.first * J] =
          -(a_ptr[0] + delta_tau * (a_ptr[1] + delta_tau * (a_ptr[2] + delta_tau * a_ptr[3])));
    }
  }

  if (I >= J_min && I < J_max)  // I==J
  {
    int t0_ind = Nt / 2;
    int r0_ind = r0_min_r1[0];

    int col_ind = b[I] + Nb * (b[I] + Nb * r0_ind);
    int row_ind = 4 * t0_ind;
    int alpha_LD = alpha_gs.first;

    assert(row_ind > -1 and row_ind < alpha_cs.first);
    assert(col_ind > -1 and col_ind < alpha_cs.second);

    assert(I > -1 and I < G0_cs.first);

    G0[I + G0_gs.first * I] = -alpha[row_ind + col_ind * alpha_LD];
  }
}

template <typename Real>
__global__ void akima_interpolation_fat_row(int Nb, int Nr, int Nt, Real beta, int Nc, int Nv,
                                            int* b, int* r, Real* t, Real* G0,
                                            std::pair<int, int> G0_cs, std::pair<int, int> G0_gs,
                                            Real* r0_min_r1, std::pair<int, int> r0_min_r1_cs,
                                            std::pair<int, int> r0_min_r1_gs, Real* alpha,
                                            std::pair<int, int> alpha_cs,
                                            std::pair<int, int> alpha_gs) {
  int I = Nc + threadIdx.x + BLOCK_SIZE_x * blockIdx.x;

  int J_min = BLOCK_SIZE_y * (blockIdx.y + 0);
  int J_max = BLOCK_SIZE_y * (blockIdx.y + 1);

  J_min = max(J_min, 0);
  J_max = min(J_max, Nc);

  if (I >= Nc && I < Nv) {
    for (int J = J_min; J < J_max; ++J) {
      int delta_r = r0_min_r1[r[J] + r[I] * r0_min_r1_gs.first];
      Real tau = t[I] - t[J];

      Real scaled_tau = (tau + beta) * Real(Nt) / (2. * beta);

      int t_ind = scaled_tau;
      Real delta_tau = scaled_tau - t_ind;

      assert(delta_tau > -1.e-16 and delta_tau < 1 + -1.e-16);

      int col_ind = b[I] + Nb * (b[J] + Nb * delta_r);
      int row_ind = 4 * t_ind;
      int alpha_LD = alpha_gs.first;

      assert(row_ind > -1 and row_ind < alpha_cs.first);
      assert(col_ind > -1 and col_ind < alpha_cs.second);

      Real* a_ptr = &alpha[row_ind + col_ind * alpha_LD];

      assert(I > -1 and I < G0_cs.first);
      assert(J > -1 and J < G0_cs.second);

      G0[I + G0_gs.first * J] =
          -(a_ptr[0] + delta_tau * (a_ptr[1] + delta_tau * (a_ptr[2] + delta_tau * a_ptr[3])));
    }
  }
}

template <class Real>
void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, Real beta, int Nc, int Nv, int* b, int* r,
                                Real* t, Real* G0, std::pair<int, int> G0_cs,
                                std::pair<int, int> G0_gs, Real* r0_min_r1,
                                std::pair<int, int> r0_min_r1_cs, std::pair<int, int> r0_min_r1_gs,
                                Real* alpha, std::pair<int, int> alpha_cs,
                                std::pair<int, int> alpha_gs) {
  // assert(cuda_check_for_errors("init 2 interpolation_kernel"));

  if (Nv - Nc > 0 and Nv > 0) {
    checkErrorsCudaDebug();

    int bl_x = dca::util::ceilDiv(Nv, BLOCK_SIZE_x);
    int bl_y = dca::util::ceilDiv(Nv - Nc, BLOCK_SIZE_y);

    dim3 threads(BLOCK_SIZE_x);
    dim3 blocks(bl_x, bl_y);

    akima_interpolation_fat_column<<<blocks, threads>>>(Nb, Nr, Nt, beta, Nc, Nv, b, r, t, G0,
                                                        G0_cs, G0_gs, r0_min_r1, r0_min_r1_cs,
                                                        r0_min_r1_gs, alpha, alpha_cs, alpha_gs);
    checkErrorsCudaDebug();
  }

  if (Nv - Nc > 0 and Nc > 0) {
    checkErrorsCudaDebug();

    int bl_x = dca::util::ceilDiv(Nv - Nc, BLOCK_SIZE_x);
    int bl_y = dca::util::ceilDiv(Nc, BLOCK_SIZE_y);

    dim3 threads(BLOCK_SIZE_x);
    dim3 blocks(bl_x, bl_y);

    akima_interpolation_fat_row<<<blocks, threads>>>(Nb, Nr, Nt, beta, Nc, Nv, b, r, t, G0, G0_cs,
                                                     G0_gs, r0_min_r1, r0_min_r1_cs, r0_min_r1_gs,
                                                     alpha, alpha_cs, alpha_gs);

    checkErrorsCudaDebug();
  }
}
template void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, float beta, int Nc, int Nv, int* b,
                                         int* r, float* t, float* G0, std::pair<int, int> G0_cs,
                                         std::pair<int, int> G0_gs, float* r0_min_r1,
                                         std::pair<int, int> r0_min_r1_cs,
                                         std::pair<int, int> r0_min_r1_gs, float* alpha,
                                         std::pair<int, int> alpha_cs, std::pair<int, int> alpha_gs);
template void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, double beta, int Nc, int Nv,
                                         int* b, int* r, double* t, double* G0,
                                         std::pair<int, int> G0_cs, std::pair<int, int> G0_gs,
                                         double* r0_min_r1, std::pair<int, int> r0_min_r1_cs,
                                         std::pair<int, int> r0_min_r1_gs, double* alpha,
                                         std::pair<int, int> alpha_cs, std::pair<int, int> alpha_gs);

template <typename Real>
void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, Real beta, int Nc, int Nv, int* b, int* r,
                                Real* t, Real* G0, std::pair<int, int> G0_cs,
                                std::pair<int, int> G0_gs, Real* r0_min_r1,
                                std::pair<int, int> r0_min_r1_cs, std::pair<int, int> r0_min_r1_gs,
                                Real* alpha, std::pair<int, int> alpha_cs,
                                std::pair<int, int> alpha_gs, int thread_id, int stream_id) {
  // assert(cuda_check_for_errors("init 2 interpolation_kernel"));

  if (Nv - Nc > 0 and Nv > 0) {
    checkErrorsCudaDebug();

    int bl_x = dca::util::ceilDiv(Nv, BLOCK_SIZE_x);
    int bl_y = dca::util::ceilDiv(Nv - Nc, BLOCK_SIZE_y);

    dim3 threads(BLOCK_SIZE_x);
    dim3 blocks(bl_x, bl_y);

    cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

    akima_interpolation_fat_column<<<blocks, threads, 0, stream_handle>>>(
        Nb, Nr, Nt, beta, Nc, Nv, b, r, t, G0, G0_cs, G0_gs, r0_min_r1, r0_min_r1_cs, r0_min_r1_gs,
        alpha, alpha_cs, alpha_gs);

    checkErrorsCudaDebug();
  }

  if (Nv - Nc > 0 and Nc > 0) {
    checkErrorsCudaDebug();

    int bl_x = dca::util::ceilDiv(Nv - Nc, BLOCK_SIZE_x);
    int bl_y = dca::util::ceilDiv(Nc, BLOCK_SIZE_y);

    dim3 threads(BLOCK_SIZE_x);
    dim3 blocks(bl_x, bl_y);

    cudaStream_t stream_handle = dca::linalg::util::getStream(thread_id, stream_id);

    akima_interpolation_fat_row<<<blocks, threads, 0, stream_handle>>>(
        Nb, Nr, Nt, beta, Nc, Nv, b, r, t, G0, G0_cs, G0_gs, r0_min_r1, r0_min_r1_cs, r0_min_r1_gs,
        alpha, alpha_cs, alpha_gs);

    checkErrorsCudaDebug();
  }
}
template void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, float beta, int Nc, int Nv, int* b,
                                         int* r, float* t, float* G0, std::pair<int, int> G0_cs,
                                         std::pair<int, int> G0_gs, float* r0_min_r1,
                                         std::pair<int, int> r0_min_r1_cs,
                                         std::pair<int, int> r0_min_r1_gs, float* alpha,
                                         std::pair<int, int> alpha_cs, std::pair<int, int> alpha_gs,
                                         int thread_id, int stream_id);
template void akima_interpolation_on_GPU(int Nb, int Nr, int Nt, double beta, int Nc, int Nv,
                                         int* b, int* r, double* t, double* G0,
                                         std::pair<int, int> G0_cs, std::pair<int, int> G0_gs,
                                         double* r0_min_r1, std::pair<int, int> r0_min_r1_cs,
                                         std::pair<int, int> r0_min_r1_gs, double* alpha,
                                         std::pair<int, int> alpha_cs, std::pair<int, int> alpha_gs,
                                         int thread_id, int stream_id);

}  // namespace g0kernels
}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca
