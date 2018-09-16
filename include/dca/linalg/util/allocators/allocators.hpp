// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides include all types of allocators, and provides a default selector.

#ifndef DCA_LINALG_UTIL_ALLOCATORS_HPP
#define DCA_LINALG_UTIL_ALLOCATORS_HPP

#include "aligned_allocator.hpp"
#include "dca/linalg/device_type.hpp"
#ifdef DCA_HAVE_CUDA
#include "device_allocator.hpp"
#include "managed_allocator.hpp"
#include "pinned_allocator.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace linalg {
namespace util {
namespace selector {
// dca::linalg::util::selector::
template <typename T, DeviceType device>
struct DefaultAllocator;

#ifdef DCA_HAVE_CUDA
template <typename T>
struct DefaultAllocator<T, CPU> {
  using type = PinnedAllocator<T>;
};

template <typename T>
struct DefaultAllocator<T, GPU> {
  using type = DeviceAllocator<T>;
};
#else

template <typename T>
struct DefaultAllocator<T, CPU> {
  using type = AlignedAllocator<T>;
};

#endif  // DCA_HAVE_CUDA

}  // selector
   // dca::linalg::util:

template <typename T, DeviceType device>
using DefaultAllocator = typename selector::DefaultAllocator<T, device>::type;

}  // util
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_ALLOCATORS_HPP
