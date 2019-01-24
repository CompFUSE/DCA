// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides typedefs for standard vectors with different allocators.

#ifndef DCA_LINALG_UTIL_VECTORS_TYPEDEFS_HPP
#define DCA_LINALG_UTIL_VECTORS_TYPEDEFS_HPP

#include <vector>

#ifdef DCA_HAVE_CUDA
#include "managed_allocator.hpp"
#include "pinned_allocator.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace linalg {
namespace util {

#ifdef DCA_HAVE_CUDA
template <typename T>
using HostVector = std::vector<T, PinnedAllocator<T>>;
template <typename T>
using ManagedVector = std::vector<T, ManagedAllocator<T>>;

#else
template <typename T>
using HostVector = std::vector<T>;
#endif  // DCA_HAVE_CUDA

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_VECTORS_TYPEDEFS_HPP
