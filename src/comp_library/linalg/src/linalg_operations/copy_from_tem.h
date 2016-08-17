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

#ifndef COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_COPY_FROM_TEM_H
#define COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_COPY_FROM_TEM_H

#include <utility>
#include "comp_library/linalg/linalg_device_types.h"

namespace LIN_ALG {

inline size_t sizeFromPair(std::pair<int,int> size) {
  return static_cast<size_t>(size.first) * static_cast<size_t>(size.second);
}

template <device_type source_device_name, device_type target_device_name>
class COPY_FROM {
  template <typename cpu_matrix_type, typename gpu_matrix_type>
  static void execute(cpu_matrix_type& cpu_matrix, gpu_matrix_type& gpu_matrix);

  template <typename cpu_matrix_type, typename gpu_matrix_type>
  static void execute(cpu_matrix_type& cpu_matrix, gpu_matrix_type& gpu_matrix, int thread_id,
                      int stream_id);

  template <typename cpu_scalartype, typename gpu_scalartype>
  static void execute(gpu_scalartype* ptr, std::pair<int, int>& c_s, std::pair<int, int>& g_s,
                      cpu_scalartype* o_ptr, std::pair<int, int>& o_c_s, std::pair<int, int>& o_g_s);
};
}

#endif  // COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_COPY_FROM_TEM_H
