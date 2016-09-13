// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Structure to copy a matrix from the GPU to the CPU.

#ifndef COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_COPY_FROM_GPU_CPU_H
#define COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_COPY_FROM_GPU_CPU_H

#include <cassert>
#include <utility>

#include "comp_library/linalg/src/linalg_operations/copy_from_tem.h"

namespace LIN_ALG {
namespace COPY_FROM_GPU_to_CPU {
// LIN_ALG::COPY_FROM_GPU_to_CPU::

template <typename scalartype>
void memcopy_d_to_h(scalartype* target_ptr, scalartype* source_ptr, int size);

template <typename scalartype>
void memcopy_d_to_h_async(scalartype* target_ptr, scalartype* source_ptr, int size, int thread_id,
                          int stream_id);

template <typename scalartype>
void memcopy_2D_d_to_h(scalartype* source_ptr, const std::pair<int, int>& source_c_s,
                       const std::pair<int, int>& source_g_s, scalartype* target_ptr,
                       const std::pair<int, int>& target_c_s, const std::pair<int, int>& target_g_s);

template <typename scalartype>
void memcopy_2D_d_to_h_async(scalartype* source_ptr, const std::pair<int, int>& source_c_s,
                             const std::pair<int, int>& source_g_s, scalartype* target_ptr,
                             const std::pair<int, int>& target_c_s, const std::pair<int, int>& target_g_s,
                             int thread_id, int stream_id);

}  // COPY_FROM_GPU_to_CPU

template <>
class COPY_FROM<GPU, CPU> {
public:
  template <typename scalartype>
  static void execute(scalartype* ptr_gpu, scalartype* ptr_cpu, int size) {
    COPY_FROM_GPU_to_CPU::memcopy_d_to_h(ptr_cpu, ptr_gpu, size);
  }

  template <typename scalartype>
  static void execute(scalartype* ptr_gpu, scalartype* ptr_cpu, int size, int thread_id,
                      int stream_id) {
    COPY_FROM_GPU_to_CPU::memcopy_d_to_h_async(ptr_cpu, ptr_gpu, size, thread_id, stream_id);
  }

  template <typename scalartype>
  static void execute(scalartype* source_ptr, const std::pair<int, int>& source_c_s,
                      const std::pair<int, int>& source_g_s, scalartype* target_ptr,
                      const std::pair<int, int>& target_c_s, const std::pair<int, int>& target_g_s) {
    assert(source_c_s == target_c_s);

    COPY_FROM_GPU_to_CPU::memcopy_2D_d_to_h(source_ptr, source_c_s, source_g_s, target_ptr,
                                            target_c_s, target_g_s);
  }

  template <typename scalartype>
  static void execute(scalartype* source_ptr, const std::pair<int, int>& source_c_s,
                      const std::pair<int, int>& source_g_s, scalartype* target_ptr,
                      const std::pair<int, int>& target_c_s, const std::pair<int, int>& target_g_s,
                      int thread_id, int stream_id) {
    assert(source_c_s == target_c_s);

    COPY_FROM_GPU_to_CPU::memcopy_2D_d_to_h_async(source_ptr, source_c_s, source_g_s, target_ptr,
                                                  target_c_s, target_g_s, thread_id, stream_id);
  }

  template <typename gpu_matrix_type, typename cpu_matrix_type>
  static void execute(gpu_matrix_type& gpu_matrix, cpu_matrix_type& cpu_matrix) {
    if (gpu_matrix.capacity() == cpu_matrix.capacity()) {
      size_t size = sizeFromPair(gpu_matrix.capacity());

      execute(gpu_matrix.ptr(), cpu_matrix.ptr(), size);
    }
    else
      execute(gpu_matrix.ptr(), gpu_matrix.size(), gpu_matrix.capacity(),
              cpu_matrix.ptr(), cpu_matrix.size(), cpu_matrix.capacity());
  }

  template <typename gpu_matrix_type, typename cpu_matrix_type>
  static void execute(gpu_matrix_type& gpu_matrix, cpu_matrix_type& cpu_matrix, int thread_id,
                      int stream_id) {
    if (gpu_matrix.capacity() == cpu_matrix.capacity()) {
      size_t size = sizeFromPair(gpu_matrix.capacity());

      execute(gpu_matrix.ptr(), cpu_matrix.ptr(), size, thread_id, stream_id);
    }
    else
      execute(gpu_matrix.ptr(), gpu_matrix.size(), gpu_matrix.capacity(),
              cpu_matrix.ptr(), cpu_matrix.size(), cpu_matrix.capacity(),
              thread_id, stream_id);
  }
};

}  // LIN_ALG

#endif  // COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_COPY_FROM_GPU_CPU_H
