// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Declares a vector that tries to allocate first on global device memory, then on managed memory.

#ifndef DCA_HAVE_CUDA
#error "This file requires CUDA."
#endif

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_VECTOR_MANAGED_FALLBACK_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_VECTOR_MANAGED_FALLBACK_HPP

#include "dca/linalg/util/allocators/device_allocator.hpp"
#include "dca/linalg/util/allocators/managed_allocator.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

template <class T>
class ManagedFallbackAllocator : public std::allocator<T> {
public:
  bool managed() const {
    return managed_;
  }

  T* allocate(std::size_t n) {
    try {
      return linalg::util::DeviceAllocator<T>().allocate(n);
    }
    catch (std::bad_alloc&) {
      managed_ = true;
      return linalg::util::ManagedAllocator<T>().allocate(n);
    }
  }

  void deallocate(T*& ptr) noexcept {
    if (managed_)
      linalg::util::ManagedAllocator<T>().deallocate(ptr);
    else
      linalg::util::DeviceAllocator<T>().deallocate(ptr);
  }

private:
  bool managed_ = false;
};

template <class T>
class VectorManagedFallback final
    : public linalg::Vector<T, linalg::GPU, ManagedFallbackAllocator<T>> {
public:
  using BaseClass = linalg::Vector<T, linalg::GPU, ManagedFallbackAllocator<T>>;

  void setStream(cudaStream_t stream) const {
    if (BaseClass::allocator_.managed())
      cudaStreamAttachMemAsync(stream, BaseClass::data());
  }
};

}  // accumulator
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_VECTOR_MANAGED_FALLBACK_HPP
