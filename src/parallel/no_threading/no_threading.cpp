// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file implements no_threading.hpp.

#include "dca/parallel/no_threading/no_threading.hpp"

namespace dca {
namespace parallel {

constexpr char NoThreading::parallel_type_str_[];

std::ostream& operator<<(std::ostream& o, const NoThreading& c) {
  o << '\n' << "threading type:" << c.parallel_type_str_ << '\n' << "number of threads:" << 1;
  return o;
}

}  // parallel
}  // dca
