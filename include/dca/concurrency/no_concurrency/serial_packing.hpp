// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides a (fake) interface to do (un)packing in serial execution.

#ifndef DCA_PARALLEL_NO_CONCURRENCY_SERIAL_PACKING_HPP
#define DCA_PARALLEL_NO_CONCURRENCY_SERIAL_PACKING_HPP

#include <cstdlib>
#include <string>
#include <vector>
#include "comp_library/function_library/function.h"

namespace dca {
namespace concurrency {
// dca::concurrency::

class SerialPacking {
public:
  template <typename T>
  std::size_t get_buffer_size(const T&) const {
    return sizeof(T);
  }
  std::size_t get_buffer_size(const std::string& str) const {
    return str.size() * sizeof(char);
  }
  template <typename T>
  std::size_t get_buffer_size(const std::vector<T>& v) const {
    return v.size() * sizeof(v[0]);
  }
  // TODO: Const correctness.
  template <typename T, class Domain>
  std::size_t get_buffer_size(FUNC_LIB::function<T, Domain>& f) const {
    return f.size() * sizeof(f(0));
  }
};

}  // concurrency
}  // dca

#endif  // DCA_PARALLEL_NO_CONCURRENCY_SERIAL_PACKING_HPP
