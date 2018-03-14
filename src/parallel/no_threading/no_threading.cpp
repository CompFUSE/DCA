// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file implements no_threading.hpp.

#include "dca/parallel/no_threading/no_threading.hpp"

namespace dca {
namespace parallel {

constexpr char NoThreading::parallel_type_str_[];

std::ostream& operator<<(std::ostream& o, const NoThreading& c) {
  o << '\n'
    << "threading type:" << c.parallel_type_str_ << '\n'
    << "number of threads:" << c.data_.num_threads;
  return o;
}

}  // parallel
}  // dca
