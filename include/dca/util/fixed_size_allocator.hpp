// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Non thread-safe allocator to quickly allocate and deallocate objects of fixes memory size.
// Inspired by https://codereview.stackexchange.com/questions/82869/fixed-size-block-allocator

#ifndef DCA_UTIL_FIXED_SIZE_ALLOCATOR_HPP
#define DCA_UTIL_FIXED_SIZE_ALLOCATOR_HPP

#include <memory>
#include <vector>

namespace dca {
namespace util {
namespace details {

template <class T, std::size_t objects_per_pool = 128>
class FixedSizeAllocator {
public:
  FixedSizeAllocator() = default;
  FixedSizeAllocator(const FixedSizeAllocator&) = delete;
  FixedSizeAllocator& operator=(const FixedSizeAllocator&) = delete;
  FixedSizeAllocator(FixedSizeAllocator&&) = default;
  FixedSizeAllocator& operator=(FixedSizeAllocator&&) = default;

  // Performs memory allocation and calls the constructor with arguments args.
  // Must be matched by a call to destroy on the pointer returned by this function.
  template <class... Args>
  [[nodiscard]] T* create(Args&&... args);

  // Calls the destructor and deallocate ptr.
  void destroy(T* ptr) noexcept;

private:
  [[nodiscard]] T* allocate();
  void deallocate(T* ptr) noexcept;
  void allocatePool();

  union TNode {
    char data[sizeof(T)];
    TNode* next;  //  pointer to the next free node.
  };

  using Pool = std::array<TNode, objects_per_pool>;

  TNode* free_ = nullptr;                     // the topmost free chunk of memory.
  std::vector<std::unique_ptr<Pool>> pools_;  // all allocated pools of memory
};

template <class T, std::size_t objects_per_pool>
template <class... Args>
T* FixedSizeAllocator<T, objects_per_pool>::create(Args&&... args) {
  T* allocation = allocate();
  return new (allocation) T(std::forward<Args>(args)...);
}

template <class T, std::size_t objects_per_pool>
void FixedSizeAllocator<T, objects_per_pool>::destroy(T* ptr) noexcept {
  if (ptr) {
    ptr->~T();
    deallocate(ptr);
  }
}

template <class T, std::size_t objects_per_pool>
T* FixedSizeAllocator<T, objects_per_pool>::allocate() {
  if (!free_) {
    allocatePool();
  }
  TNode* result = free_;  // allocate the topmost element.
  free_ = free_->next;    // and pop it from the stack of free chunks
  return reinterpret_cast<T*>(result);
}

template <class T, std::size_t objects_per_pool>
void FixedSizeAllocator<T, objects_per_pool>::deallocate(T* ptr) noexcept {
  TNode* node = reinterpret_cast<TNode*>(ptr);
  // add to the stack of chunks
  node->next = free_;
  free_ = node;
}

template <class T, std::size_t objects_per_pool>
void FixedSizeAllocator<T, objects_per_pool>::allocatePool() {
  // Allocate new memory.
  pools_.emplace_back(std::make_unique<Pool>());

  // Form a stack from this pool.
  auto& new_pool = *pools_.back();
  for (int i = 0; i < objects_per_pool - 1; ++i) {
    new_pool[i].next = &new_pool[i + 1];
  }
  new_pool.back().next = nullptr;

  free_ = new_pool.data();
}

}  // namespace details
}  // namespace util
}  // namespace dca

#endif  // DCA_UTIL_FIXED_SIZE_ALLOCATOR_HPP
