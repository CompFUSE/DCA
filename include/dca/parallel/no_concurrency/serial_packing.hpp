// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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

#include "dca/function/function.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

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
  std::size_t get_buffer_size(func::function<T, Domain>& f) const {
    return f.size() * sizeof(f(0));
  }
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_NO_CONCURRENCY_SERIAL_PACKING_HPP
