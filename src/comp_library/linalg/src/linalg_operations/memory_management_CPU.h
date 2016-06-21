// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_MEMORY_MANAGEMENT_CPU_H
#define COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_MEMORY_MANAGEMENT_CPU_H

#include <cassert>
#include <iostream>
#include <utility>

#include "comp_library/linalg/src/linalg_operations/memory_management_tem.h"

namespace LIN_ALG {
namespace MEMORY_MANAGEMENT_ON_GPU {
// LIN_ALG::MEMORY_MANAGEMENT_ON_GPU::

template <typename scalartype>
void allocate_pinned_host_memory(scalartype*& ptr, int global_size);
template <typename scalartype>
void allocate_pinned_host_memory(scalartype*& ptr, std::pair<int, int> global_size);

template <typename scalartype>
void deallocate_pinned_host_memory(scalartype*& ptr);

}  // MEMORY_MANAGEMENT_ON_GPU

template <>
class MEMORY_MANAGEMENT<CPU> {
public:
  template <typename scalartype>
  inline static scalartype& get(scalartype* ptr) {
    return ptr[0];
  }

  template <typename scalartype>
  inline static scalartype& get(scalartype* ptr, int index) {
    return ptr[index];
  }

  template <typename scalartype>
  inline static void set(scalartype* ptr, scalartype val) {
    ptr[0] = val;
  }

  template <typename scalartype>
  inline static void add(scalartype* ptr, scalartype val) {
    ptr[0] += val;
  }

  template <typename scalartype>
  inline static void allocate(scalartype*& ptr, int global_size) {
    assert(ptr == NULL);

#ifdef ENABLE_PINNED_MEMORY_ALLOCATION
    MEMORY_MANAGEMENT_ON_GPU::allocate_pinned_host_memory(ptr, global_size);
#else
    posix_memalign((void**)&ptr, 128, global_size * sizeof(scalartype));
#endif

    assert(ptr != NULL);
  }

  template <typename scalartype>
  inline static void allocate(scalartype*& ptr, std::pair<int, int>& global_size) {
    assert(ptr == NULL);

#ifdef ENABLE_PINNED_MEMORY_ALLOCATION
    MEMORY_MANAGEMENT_ON_GPU::allocate_pinned_host_memory(ptr, global_size);
#else
    posix_memalign((void**)&ptr, 128, global_size.first * global_size.second * sizeof(scalartype));
#endif

    assert(ptr != NULL);
  }

  template <typename scalartype>
  inline static void deallocate(scalartype*& ptr) {
    assert(ptr != NULL);

#ifdef ENABLE_PINNED_MEMORY_ALLOCATION
    MEMORY_MANAGEMENT_ON_GPU::deallocate_pinned_host_memory(ptr);
#else
    free(ptr);
#endif

    ptr = NULL;
  }

  /*
    template<typename scalartype>
    inline static void memcopy(scalartype* target_ptr, scalartype* source_ptr, int size){
        memcpy(target_ptr, source_ptr, sizeof(scalartype)*size);
    }
  */

  template <typename scalartype>
  inline static void set_to_zero(scalartype* ptr, int size) {
    for (int l = 0; l < size; ++l)
      ptr[l] = scalartype(0);
  }

  template <typename scalartype>
  inline static void set_to_zero(scalartype* ptr, int LD, int size) {
    for (int l = 0; l < size; ++l)
      ptr[l * LD] = scalartype(0);
  }

  template <typename scalartype>
  static void print(scalartype* ptr, int c_s, int g_s) {
    std::cout.precision(6);
    std::cout << std::scientific;

    std::cout << "\n\n";
    std::cout << "\t current-size : " << c_s << "\n";
    std::cout << "\t global -size : " << g_s << "\n";
    std::cout << "\n\n";

    for (int i = 0; i < c_s; i++)
      std::cout << "\t" << ptr[i];
    std::cout << "\n";

    std::cout << "\n\n\n";
  }

  template <typename scalartype>
  static void print(scalartype*& ptr, std::pair<int, int>& c_s, std::pair<int, int>& g_s) {
    std::cout.precision(6);
    std::cout << std::scientific;

    std::cout << "\n\n";
    std::cout << "\t current-size : " << c_s.first << "\t" << c_s.second << "\n";
    std::cout << "\t global -size : " << g_s.first << "\t" << g_s.second << "\n";
    std::cout << "\n\n";

    for (int i = 0; i < c_s.first; i++) {
      for (int j = 0; j < c_s.second; j++)
        std::cout << "\t" << ptr[i + g_s.first * j];
      std::cout << "\n";
    }

    std::cout << "\n\n\n";
  }
};

}  // LIN_ALG

#endif  // COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_MEMORY_MANAGEMENT_CPU_H
