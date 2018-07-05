// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Provides an utility to select block and grid size for CUDA kernels.

#ifndef DCA_UTIL_CUDA_BLOCKS_HPP
#define DCA_UTIL_CUDA_BLOCKS_HPP

#include <array>
#include <cassert>

#include <cuda.h>

#include "dca/util/integer_division.hpp"

namespace dca {
namespace util {
// dca::util::

inline std::array<int, 2> get1DBlockSize(const int ni, const int block_size) {
  const int n_threads = std::min(block_size, ni);
  const int n_blocks = util::ceilDiv(ni, n_threads);
  return std::array<int, 2>{n_blocks, n_threads};
}

inline std::array<dim3, 2> get2DBlockSize(const uint i, const uint j, const uint block_size = 32) {
  assert(i > 0 && j > 0);
  const uint n_threads_i = std::min(block_size, i);
  const uint n_threads_j = std::min(block_size, j);
  if (n_threads_i * n_threads_j > 32 * 32)
    throw(std::logic_error("Block size is too big"));

  const uint n_blocks_i = dca::util::ceilDiv(i, n_threads_i);
  const uint n_blocks_j = dca::util::ceilDiv(j, n_threads_j);

  return std::array<dim3, 2>{dim3(n_blocks_i, n_blocks_j), dim3(n_threads_i, n_threads_j)};
}

inline std::array<dim3, 2> getBlockSize3D(const uint i, const uint j, const uint k) {
  const uint n_threads_k = std::min(uint(8), k);
  const uint max_block_size_ij = n_threads_k > 1 ? 8 : 32;
  const uint n_threads_i = std::min(max_block_size_ij, i);
  const uint n_threads_j = std::min(max_block_size_ij, j);

  const uint n_blocks_i = dca::util::ceilDiv(i, n_threads_i);
  const uint n_blocks_j = dca::util::ceilDiv(j, n_threads_j);
  const uint n_blocks_k = dca::util::ceilDiv(k, n_threads_k);

  return std::array<dim3, 2>{dim3(n_blocks_i, n_blocks_j, n_blocks_k),
                             dim3(n_threads_i, n_threads_j, n_blocks_k)};
}

inline auto getBlockSize(const uint i, const uint j, const uint block_size = 32) {
  return get2DBlockSize(i, j, block_size);
}

}  // util
}  // dca

#endif  // DCA_UTIL_CUDA_BLOCKS_HPP
