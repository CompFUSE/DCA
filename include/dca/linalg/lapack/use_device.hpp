// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides device-generic wrappers for some of the LAPACK routines.

#ifndef DCA_LINALG_LAPACK_USE_DEVICE_HPP
#define DCA_LINALG_LAPACK_USE_DEVICE_HPP

#include "dca/linalg/device_type.hpp"

#include "dca/linalg/lapack/lapack.hpp"

#ifdef DCA_HAVE_CUDA
#include "dca/linalg/lapack/laset_gpu.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace linalg {
namespace lapack {
// dca::linalg::lapack::

template <DeviceType DeviceName>
struct UseDevice;

template <>
struct UseDevice<CPU> {
  // Auxiliary routines
  template <typename ScalarType>
  inline static void laset(int m, int n, ScalarType offdiag, ScalarType diag, ScalarType* a,
                           int lda, int /*thread_id*/, int /*stream_id*/) {
    lapack::laset("A", m, n, offdiag, diag, a, lda);
  }
};

#ifdef DCA_HAVE_CUDA
template <>
struct UseDevice<GPU> {
  // Auxiliary routines
  template <typename ScalarType>
  inline static void laset(int m, int n, ScalarType offdiag, ScalarType diag, ScalarType* a,
                           int lda, int thread_id, int stream_id) {
    lapack::laset_gpu(m, n, offdiag, diag, a, lda, thread_id, stream_id);
  }
};
#endif  // DCA_HAVE_CUDA

}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_USE_DEVICE_HPP
