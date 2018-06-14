// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides a cublas handle container.

#ifndef DCA_LINALG_UTIL_HANDLE_CONTAINER_HPP
#define DCA_LINALG_UTIL_HANDLE_CONTAINER_HPP

#include <array>
#include <cassert>
#include <cublas_v2.h>
#include "dca/linalg/util/error_cublas.hpp"
#include "dca/linalg/util/stream_container.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

template <size_t max_threads>
class HandleContainer {
public:
  HandleContainer() {
    for (size_t i = 0; i < handles_.size(); ++i) {
      cublasStatus_t ret = cublasCreate(&(handles_[i]));
      checkRC(ret);
    }
  }

  // Initializes the handels and sets the stream of handle(id) with stream(id, 0)
  // for 0 <= id < max_threads.
  template <size_t streams_per_thread>
  HandleContainer(StreamContainer<max_threads, streams_per_thread>& stream_container)
      : HandleContainer() {
    for (size_t i = 0; i < handles_.size(); ++i)
      setStream(i, stream_container(i, 0));
  }

  ~HandleContainer() {
    for (size_t i = 0; i < handles_.size(); ++i) {
      cublasStatus_t ret = cublasDestroy(handles_[i]);
      checkRC(ret);
    }
  }

  HandleContainer(const HandleContainer<max_threads>&) = delete;
  HandleContainer& operator=(const HandleContainer<max_threads>&) = delete;

  // Returns the cublasHandle associated with thread 'thread_id'.
  // Preconditions: 0 <= thread_id < DCA_MAX_THREADS.
  cublasHandle_t operator()(int thread_id) {
    assert(thread_id >= 0 && static_cast<size_t>(thread_id) < max_threads);
    return handles_[thread_id];
  }

  // Sets the stream of the handle associated with thread 'thread_id'.
  // Preconditions: 0 <= thread_id < DCA_MAX_THREADS.
  void setStream(int thread_id, cudaStream_t stream) {
    cublasStatus_t ret = cublasSetStream(handles_[thread_id], stream);
    checkRC(ret);
  }

private:
  std::array<cublasHandle_t, max_threads> handles_;
};
}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_HANDLE_CONTAINER_HPP
