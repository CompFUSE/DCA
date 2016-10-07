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

#include "comp_library/linalg/src/linalg_operations/memory_management_tem.h"

namespace LIN_ALG {
namespace MEMORY_MANAGEMENT_ON_GPU {
// LIN_ALG::MEMORY_MANAGEMENT_ON_GPU::

template <typename scalartype>
void remove_first_row(int m, int n, scalartype* A, int LDA);
template <typename scalartype>
void remove_first_col(int m, int n, scalartype* A, int LDA);

}  // MEMORY_MANAGEMENT_ON_GPU

template <>
class MEMORY_MANAGEMENT<GPU> {
public:

  template <typename scalartype>
  static void remove_first_row(int m, int n, scalartype* A, int LDA) {
    MEMORY_MANAGEMENT_ON_GPU::remove_first_row(m, n, A, LDA);
  }

  template <typename scalartype>
  static void remove_first_col(int m, int n, scalartype* A, int LDA) {
    MEMORY_MANAGEMENT_ON_GPU::remove_first_col(m, n, A, LDA);
  }

};
}

#endif  // COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_MEMORY_MANAGEMENT_GPU_H
