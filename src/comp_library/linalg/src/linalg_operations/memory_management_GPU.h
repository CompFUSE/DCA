// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_MEMORY_MANAGEMENT_GPU_H
#define COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_MEMORY_MANAGEMENT_GPU_H

#include <iostream>
#include <utility>

#include "comp_library/linalg/src/linalg_operations/copy_from_GPU_CPU.h"
#include "comp_library/linalg/src/linalg_operations/memory_management_tem.h"

namespace LIN_ALG {
namespace MEMORY_MANAGEMENT_ON_GPU {
// LIN_ALG::MEMORY_MANAGEMENT_ON_GPU::

template <typename scalartype>
scalartype get(scalartype* ptr);
template <typename scalartype>
scalartype get(scalartype* ptr, int index);

template <typename scalartype>
void set(scalartype* ptr, scalartype val);
template <typename scalartype>
void add(scalartype* ptr, scalartype val);

template <typename scalartype>
void allocate(scalartype*& ptr, int global_size);
template <typename scalartype>
void allocate(scalartype*& ptr, std::pair<int, int> global_size);
template <typename scalartype>
void deallocate(scalartype*& ptr);

//     template<typename scalartype> void memcopy_d_to_d(scalartype* target_ptr, scalartype*
//     source_ptr, int size);
//     template<typename scalartype> void memcopy_d_to_h(scalartype* target_ptr, scalartype*
//     source_ptr, int size);
//     template<typename scalartype> void memcopy_h_to_d(scalartype* target_ptr, scalartype*
//     source_ptr, int size);

//     template<typename scalartype> void memcopy_d_to_d_async(scalartype* target_ptr, scalartype*
//     source_ptr, int size);
//     template<typename scalartype> void memcopy_d_to_h_async(scalartype* target_ptr, scalartype*
//     source_ptr, int size);
//     template<typename scalartype> void memcopy_h_to_d_async(scalartype* target_ptr, scalartype*
//     source_ptr, int size);

template <typename scalartype>
void set_to_zero(scalartype* ptr, int m);
template <typename scalartype>
void set_to_zero(scalartype* ptr, int LD, int m);

template <typename scalartype>
void remove_first_row(int m, int n, scalartype* A, int LDA);
template <typename scalartype>
void remove_first_col(int m, int n, scalartype* A, int LDA);

}  // MEMORY_MANAGEMENT_ON_GPU

template <>
class MEMORY_MANAGEMENT<GPU> {
public:
  template <typename scalartype>
  inline static scalartype get(scalartype* ptr) {
    return MEMORY_MANAGEMENT_ON_GPU::get(ptr);
  }

  template <typename scalartype>
  inline static scalartype get(scalartype* ptr, int index) {
    return MEMORY_MANAGEMENT_ON_GPU::get(ptr, index);
  }

  template <typename scalartype>
  inline static void set(scalartype* ptr, scalartype val) {
    MEMORY_MANAGEMENT_ON_GPU::set(ptr, val);
  }

  template <typename scalartype>
  inline static void add(scalartype* ptr, scalartype val) {
    MEMORY_MANAGEMENT_ON_GPU::add(ptr, val);
  }

  template <typename scalartype>
  static void allocate(scalartype*& ptr, int global_size) {
    MEMORY_MANAGEMENT_ON_GPU::allocate(ptr, global_size);
  }

  template <typename scalartype>
  static void allocate(scalartype*& ptr, std::pair<int, int> global_size) {
    MEMORY_MANAGEMENT_ON_GPU::allocate(ptr, global_size);
  }

  template <typename scalartype>
  static void deallocate(scalartype*& ptr) {
    MEMORY_MANAGEMENT_ON_GPU::deallocate(ptr);
  }

  template <typename scalartype>
  static void set_to_zero(scalartype* ptr, int size) {
    MEMORY_MANAGEMENT_ON_GPU::set_to_zero(ptr, size);
  }

  template <typename scalartype>
  static void set_to_zero(scalartype* ptr, int LD, int size) {
    MEMORY_MANAGEMENT_ON_GPU::set_to_zero(ptr, LD, size);
  }

  template <typename scalartype>
  static void remove_first_row(int m, int n, scalartype* A, int LDA) {
    MEMORY_MANAGEMENT_ON_GPU::remove_first_row(m, n, A, LDA);
  }

  template <typename scalartype>
  static void remove_first_col(int m, int n, scalartype* A, int LDA) {
    MEMORY_MANAGEMENT_ON_GPU::remove_first_col(m, n, A, LDA);
  }

  template <typename scalartype>
  static void print(scalartype* ptr, int c_s, int g_s) {
    int SIZE = g_s;
    scalartype* new_data = new scalartype[SIZE];

    // MEMORY_MANAGEMENT_ON_GPU::memcopy_d_to_h(new_data, ptr, SIZE);
    COPY_FROM<GPU, CPU>::execute(ptr, new_data, SIZE);

    std::cout.precision(6);
    std::cout << std::scientific;

    std::cout << "\n\n";
    std::cout << "\t current-size : " << c_s << "\n";
    std::cout << "\t global -size : " << g_s << "\n";
    std::cout << "\n\n";

    for (int i = 0; i < c_s; i++)
      std::cout << "\t" << new_data[i];
    std::cout << "\n";

    std::cout << "\n\n\n";

    delete[] new_data;
  }

  template <typename scalartype>
  static void print(scalartype* ptr, std::pair<int, int>& c_s, std::pair<int, int>& g_s) {
    int SIZE = g_s.first * g_s.second;
    scalartype* new_data = new scalartype[SIZE];

    // MEMORY_MANAGEMENT_ON_GPU::memcopy_d_to_h(new_data, ptr, SIZE);
    COPY_FROM<GPU, CPU>::execute(ptr, new_data, SIZE);

    std::cout.precision(6);
    std::cout << std::scientific;

    std::cout << "\n\n";
    std::cout << "\t current-size : " << c_s.first << "\t" << c_s.second << "\n";
    std::cout << "\t global -size : " << g_s.first << "\t" << g_s.second << "\n";
    std::cout << "\n\n";

    for (int i = 0; i < c_s.first; i++) {
      for (int j = 0; j < c_s.second; j++)
        std::cout << "\t" << new_data[i + g_s.first * j];
      std::cout << "\n";
    }

    std::cout << "\n\n\n";

    delete[] new_data;
  }
};
}

#endif  // COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_MEMORY_MANAGEMENT_GPU_H
