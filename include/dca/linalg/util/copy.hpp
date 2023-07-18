// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides memory copy utilities.

#ifndef DCA_LINALG_UTIL_COPY_HPP
#define DCA_LINALG_UTIL_COPY_HPP

#include <cassert>
#include <complex>
#include <cstring>
#include "dca/linalg/device_type.hpp"
#include "gpu_stream.hpp"

#ifdef DCA_HAVE_GPU
#include "dca/platform/dca_gpu.h"
#endif
#include "dca/linalg/util/stream_functions.hpp"

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

template <typename ScalarType>
inline void memoryCopyCpu(ScalarType* dest, const ScalarType* src, size_t sz) {
  std::memcpy(dest, src, sz * sizeof(ScalarType));
}

template <typename T>
const void* constCastToVoid(T t) {
  static_assert(
      std::is_const_v<typename std::remove_pointer<typename std::remove_reference<T>::type>::type>);
  if constexpr (std::is_array_v<std::remove_pointer<T>>)
    return static_cast<const void*>(&t[0]);
  else
    return static_cast<const void*>(t);
}

template <typename T>
void* castToVoid(T t) {
  if constexpr (std::is_array_v<typename std::remove_pointer<typename std::remove_reference<T>::type>::type>)
    return (void*)t;
  else
    return static_cast<void*>(t);
}

template <typename Scalar1, typename Scalar2>
inline void memoryCopyCpu(Scalar1* dest, const Scalar2* src, size_t sz) {
  static_assert(sizeof(Scalar1) == sizeof(Scalar2), "can't copy between unequally signed types");
  std::memcpy(castToVoid(dest), constCastToVoid(src), sz * sizeof(Scalar1));
}

template <typename ScalarType>
void memoryCopyCpu(ScalarType* dest, int ld_dest, const ScalarType* src, int ld_src,
                   std::pair<int, int> size) {
  assert(size.first <= ld_dest);
  assert(size.first <= ld_src);
  assert(size.first >= 0);
  assert(size.second >= 0);
  size_t ncols = size.second;
  for (size_t i = 0; i < ncols; ++i) {
    memoryCopyCpu(dest + i * ld_dest, src + i * ld_src, size.first);
  }
}

#ifdef DCA_HAVE_GPU
// Fully synchronous 1D memory copy, i.e. all operations in the GPU queue are executed before the
// execution of this copy.
// The host continues the execution of the program when the copy is terminated.
template <typename ScalarType>
void memoryCopy(ScalarType* dest, const ScalarType* src, size_t size) {
  if (size == 0)
    return;
  cudaError_t ret = cudaMemcpy(dest, src, size * sizeof(ScalarType), cudaMemcpyDefault);
  checkRC(ret);
}

// Fully synchronous 1D memory copy, i.e. all operations in the GPU queue are executed before the
// execution of this copy.
// The host continues the execution of the program when the copy is terminated.
template <typename Scalar1, typename Scalar2>
void memoryCopy(Scalar1* dest, const Scalar2* src, size_t size) {
  if (size == 0)
    return;
  static_assert(sizeof(Scalar1) == sizeof(Scalar2));
  cudaError_t ret =
      cudaMemcpy(castToVoid(dest), constCastToVoid(src), size * sizeof(Scalar1), cudaMemcpyDefault);
  checkRC(ret);
}

template <typename Scalar1, typename Scalar2>
void memoryCopyH2H(Scalar1* dest, const Scalar2* src, size_t size) {
  if (size == 0)
    return;
  static_assert(sizeof(Scalar1) == sizeof(Scalar2));
  cudaError_t ret = cudaMemcpy(castToVoid(dest), constCastToVoid(src), size * sizeof(Scalar1),
                               cudaMemcpyHostToHost);
  checkRC(ret);
}

template <typename Scalar1, typename Scalar2>
void memoryCopyH2D(Scalar1* dest, const Scalar2* src, size_t size) {
  if (size == 0)
    return;
  static_assert(sizeof(Scalar1) == sizeof(Scalar2));
  cudaError_t ret =
      cudaMemcpy(castToVoid(dest), castToVoid(src), size * sizeof(Scalar1), cudaMemcpyHostToDevice);
  checkRC(ret);
}

// Fully synchronous 1D memory copy, i.e. all operations in the GPU queue are executed before the
// execution of this copy.
// The host continues the execution of the program when the copy is terminated.
template <typename Scalar1, typename Scalar2>
void memoryCopyD2H(Scalar1* dest, const Scalar2* src, size_t size) {
  if (size == 0)
    return;
  static_assert(sizeof(Scalar1) == sizeof(Scalar2));
  cudaError_t ret = cudaMemcpy(castToVoid(dest), constCastToVoid(src), size * sizeof(Scalar1),
                               cudaMemcpyDeviceToHost);
  checkRC(ret);
}

// Fully synchronous 2D memory copy, i.e. all operations in the GPU queue are executed before the
// execution of this copy.
// The host continues the execution of the program when the copy is terminated.
// Preconditions: ld_dest >= size.first, ld_src >= size.first.
template <typename ScalarType>
void memoryCopy(ScalarType* dest, int ld_dest, const ScalarType* src, int ld_src,
                std::pair<int, int> size) {
  if (ld_dest == 0 || ld_src == 0 || (size.first == 0 && size.second == 0))
    return;
  cudaError_t ret = cudaMemcpy2D(dest, ld_dest * sizeof(ScalarType), src, ld_src * sizeof(ScalarType),
                                 size.first * sizeof(ScalarType), size.second, cudaMemcpyDefault);
  try {
    checkRC(ret);
  }
  catch (...) {
    std::cout << "Failed memorycopy!\n";
    throw;
  }
}

template <typename Scalar1, typename Scalar2>
void memoryCopy(Scalar1* dest, int ld_dest, const Scalar2* src, int ld_src, std::pair<int, int> size) {
  static_assert(sizeof(Scalar1) == sizeof(Scalar2));
  if (ld_dest == 0 || ld_src == 0 || (size.first == 0 && size.second == 0))
    return;
  cudaError_t ret = cudaMemcpy2D(castToVoid(dest), ld_dest * sizeof(Scalar1), constCastToVoid(src),
                                 ld_src * sizeof(Scalar2), size.first * sizeof(Scalar1),
                                 size.second, cudaMemcpyDefault);
  try {
    checkRC(ret);
  }
  catch (...) {
    std::cout << "Failed memorycopy!\n";
    throw;
  }
}

template <typename Scalar1, typename Scalar2>
void memoryCopyH2D(Scalar1* dest, int ld_dest, const Scalar2* src, int ld_src,
                   std::pair<int, int> size) {
  static_assert(sizeof(Scalar1) == sizeof(Scalar2));
  if (ld_dest == 0 || ld_src == 0 || (size.first == 0 && size.second == 0))
    return;
  cudaError_t ret = cudaMemcpy2D(castToVoid(dest), ld_dest * sizeof(Scalar1), constCastToVoid(src),
                                 ld_src * sizeof(Scalar2), size.first * sizeof(Scalar1),
                                 size.second, cudaMemcpyHostToDevice);
  try {
    checkRC(ret);
  }
  catch (...) {
    std::cout << "Failed memorycopy!\n";
    throw;
  }
}

template <typename Scalar1, typename Scalar2>
void memoryCopyD2H(Scalar1* dest, int ld_dest, const Scalar2* src, int ld_src,
                   std::pair<int, int> size) {
  static_assert(sizeof(Scalar1) == sizeof(Scalar2));
  if (ld_dest == 0 || ld_src == 0 || (size.first == 0 && size.second == 0))
    return;
  cudaError_t ret = cudaMemcpy2D(castToVoid(dest), ld_dest * sizeof(Scalar1), constCastToVoid(src),
                                 ld_src * sizeof(Scalar2), size.first * sizeof(Scalar1),
                                 size.second, cudaMemcpyDeviceToHost);
  try {
    checkRC(ret);
  }
  catch (...) {
    std::cout << "Failed memorycopy!\n";
    throw;
  }
}

// Asynchronous 1D memory copy.
template <typename ScalarType>
void memoryCopyAsync(ScalarType* dest, const ScalarType* src, size_t size, const cudaStream_t stream) {
  if (size == 0)
    return;
  cudaError_t ret = cudaMemcpyAsync(dest, src, size * sizeof(ScalarType), cudaMemcpyDefault, stream);
  try {
    checkRC(ret);
  }
  catch (...) {
    std::cout << "Failed memorycopy!\n";
    throw;
  }
}

// Asynchronous 1D memory copy (stream = getStream(thread_id, stream_id)).
// Preconditions: 0 <= thread_id < DCA_MAX_THREADS,
//                0 <= stream_id < DCA_STREAMS_PER_THREADS.
template <typename ScalarType>
void memoryCopyAsync(ScalarType* dest, const ScalarType* src, size_t size, int thread_id,
                     int stream_id = 0) {
  memoryCopyAsync(dest, src, size, getStream(thread_id, stream_id));
}

// Asynchronous 2D memory copy.
// Preconditions: ld_dest >= size.first, ld_src >= size.firs.
template <typename ScalarType>
void memoryCopyAsync(ScalarType* dest, int ld_dest, const ScalarType* src, int ld_src,
                     std::pair<int, int> size, const cudaStream_t stream) {
  if (ld_dest == 0 || ld_src == 0 || (size.first == 0 && size.second == 0))
    return;
  cudaError_t ret =
      cudaMemcpy2DAsync(dest, ld_dest * sizeof(ScalarType), src, ld_src * sizeof(ScalarType),
                        size.first * sizeof(ScalarType), size.second, cudaMemcpyDefault, stream);
  try {
    checkRC(ret);
  }
  catch (...) {
    std::cout << "Failed memorycopy!\n";
    throw;
  }
}

// Asynchronous 2D memory copy (stream = getStream(thread_id, stream_id)).
// Preconditions: ld_dest >= size.first, ld_src >= size.first,
//                0 <= thread_id < DCA_MAX_THREADS,
//                0 <= stream_id < DCA_STREAMS_PER_THREADS.
template <typename ScalarType>
void memoryCopyAsync(ScalarType* dest, int ld_dest, const ScalarType* src, int ld_src,
                     std::pair<int, int> size, int thread_id, int stream_id) {
  memoryCopyAsync(dest, ld_dest, src, ld_src, size, getStream(thread_id, stream_id));
}

// Asynchronous 1D memory copy (stream = getStream(thread_id, stream_id))
// + synchronization of stream.
// Preconditions: 0 <= thread_id < DCA_MAX_THREADS,
template <typename ScalarType>
void memoryCopy(ScalarType* dest, const ScalarType* src, size_t size, int thread_id, int stream_id) {
  memoryCopyAsync(dest, src, size, thread_id, stream_id);
  syncStream(thread_id, stream_id);
}

// Asynchronous 2D memory copy (stream = getStream(thread_id, stream_id))
// + synchronization of stream.
// Preconditions: ld_dest >= size.first, ld_src >= size.first,
//                0 <= thread_id < DCA_MAX_THREADS,
//                0 <= stream_id < DCA_STREAMS_PER_THREADS.
template <typename ScalarType>
void memoryCopy(ScalarType* dest, int ld_dest, const ScalarType* src, int ld_src,
                std::pair<int, int> size, int thread_id, int stream_id) {
  memoryCopyAsync(dest, ld_dest, src, ld_src, size, thread_id, stream_id);
  syncStream(thread_id, stream_id);
}

#else

template <typename ScalarType>
void memoryCopy(ScalarType* dest, const ScalarType* src, size_t sz, int /*thread_id*/ = 0,
                int /*stream_id*/ = 0) {
  memoryCopyCpu(dest, src, sz);
}

template <typename Scalar1, typename Scalar2>
void memoryCopy(Scalar1* dest, const Scalar2* src, size_t sz, int /*thread_id*/ = 0,
                int /*stream_id*/ = 0) {
  memoryCopyCpu(dest, src, sz);
}

template <typename ScalarType>
void memoryCopy(ScalarType* dest, int ld_dest, const ScalarType* src, int ld_src,
                std::pair<int, int> size, int /*thread_id*/ = 0, int /*stream_id*/ = 0) {
  memoryCopyCpu(dest, ld_dest, src, ld_src, size);
}

template <typename Scalar1, typename Scalar2>
void memoryCopy(Scalar1* dest, int ld_dest, const Scalar2* src, int ld_src,
                std::pair<int, int> size, int /*thread_id*/ = 0, int /*stream_id*/ = 0) {
  static_assert(sizeof(Scalar1) == sizeof(Scalar2));
  memoryCopyCpu(dest, ld_dest, src, ld_src, size);
}

// Synchronous 1D memory copy fallback.
template <typename ScalarType>
void memoryCopyAsync(ScalarType* dest, const ScalarType* src, size_t size,
                     const util::GpuStream& /*s*/) {
  memoryCopyCpu(dest, src, size);
}
template <typename ScalarType>
void memoryCopyAsync(ScalarType* dest, int ld_dest, const ScalarType* src, int ld_src,
                     std::pair<int, int> size, const util::GpuStream& /*s*/) {
  memoryCopyCpu(dest, ld_dest, src, ld_src, size);
}

template <typename Scalar1, typename Scalar2>
void memoryCopyH2D(Scalar1* dest, int ld_dest, const Scalar2* src, int ld_src,
                   std::pair<int, int> size) {
  throw std::runtime_error("memoryCopyH2D should never be called in a non GPU build.");
}

template <typename Scalar1, typename Scalar2>
void memoryCopyD2H(Scalar1* dest, int ld_dest, const Scalar2* src, int ld_src,
                   std::pair<int, int> size) {
  throw std::runtime_error("memoryCopyH2D should never be called in a non GPU build.");
}

#endif  // DCA_HAVE_GPU

}  // namespace util
}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_UTIL_COPY_HPP
