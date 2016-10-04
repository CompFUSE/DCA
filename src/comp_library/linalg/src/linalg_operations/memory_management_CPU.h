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

#ifndef COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_MEMORY_MANAGEMENT_CPU_H
#define COMP_LIBRARY_LINALG_SRC_LINALG_OPERATIONS_MEMORY_MANAGEMENT_CPU_H

#include <cassert>
#include <iostream>
#include <utility>

#include "comp_library/linalg/src/linalg_operations/memory_management_tem.h"

namespace LIN_ALG {
template <>
class MEMORY_MANAGEMENT<CPU> {
  public:

  template <typename scalartype>
  static void print(const scalartype* ptr, int c_s, int g_s) {
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
  static void print(const scalartype* ptr, std::pair<int, int> c_s, std::pair<int, int> g_s) {
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
