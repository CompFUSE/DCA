// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the device methods of CtintWalker.

#include "dca/phys/dca_step/cluster_solver/ctint/walker/kernels_interface.hpp"

#include <cassert>
#include <cuda.h>
#include <cuda_runtime.h>

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {

void __global__ smallInverseKernel(const MatrixView m_inp, MatrixView m_out) {
  const int n = m_inp.nrRows();
  //  const int i = threadIdx.x;
  //  const int j = threadIdx.y;

  switch (n) {
    case 1:
      m_out(0, 0) = 1. / m_inp(0, 0);
      break;
    case 2:
      const double det = m_inp(0, 0) * m_inp(1, 1) - m_inp(1, 0) * m_inp(0, 1);
      m_out(0, 0) = m_inp(1, 1) / det;
      m_out(0, 1) = -m_inp(0, 1) / det;
      m_out(1, 0) = -m_inp(1, 0) / det;
      m_out(1, 1) = m_inp(0, 0) / det;
      break;
    default:  // abort
      asm("trap;");
  }
}

void __global__ smallInverseKernel(MatrixView m) {
  const int n = m.nrRows();
  //  const int i = threadIdx.x;
  //  const int j = threadIdx.y;
  auto swap = [](double& a, double& b) {
    const double tmp = a;
    a = b;
    b = tmp;
  };

  switch (n) {
    case 1:
      m(0, 0) = 1. / m(0, 0);
      break;
    case 2:
      const double det = m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1);
      m(0, 0) /= det;
      m(1, 0) /= -det;
      m(0, 1) /= -det;
      m(1, 1) /= det;

      swap(m(0, 0), m(1, 1));
      swap(m(1, 0), m(0, 1));
      break;
    default:  // abort
      asm("trap;");
  }
}

void smallInverse(const MatrixView& m_inp, MatrixView& m_out, cudaStream_t stream) {
  assert(m_inp.is_square() && m_out.is_square());
  const int n = m_inp.nrCols();
  smallInverseKernel<<<1, 1 /*dim3(n, n)*/, 0, stream>>>(m_inp, m_out);
}

void smallInverse(MatrixView& in_out, cudaStream_t stream) {
  assert(in_out.is_square());
  //  const int n = in_out.nrCols();
  smallInverseKernel<<<1, 1 /*dim3(n, n)*/, 0, stream>>>(in_out);
}

void __global__ separateIndexDeterminantKernel(const MatrixView M, const ushort* indices,
                                               const int n, double* det) {
  switch (n) {
    case 1:
      *det = M(indices[0], indices[0]);
      break;
    case 2:
      *det = M(indices[0], indices[0]) * M(indices[1], indices[1]) -
             M(indices[1], indices[0]) * M(indices[0], indices[1]);
      break;
    default:  // abort
      asm("trap;");
  }
}

double separateIndexDeterminant(MatrixView m, const ushort* indices, int n_indices, double* det,
                                cudaStream_t stream) {
  separateIndexDeterminantKernel<<<1, 1, 0, stream>>>(m, indices, n_indices, det);
  double host_det;
  cudaMemcpyAsync(&host_det, det, sizeof(double), cudaMemcpyDeviceToHost, stream);
  cudaStreamSynchronize(stream);
  return host_det;
}

}  // details
}  // ctint
}  // solver
}  // phys
}  // dca