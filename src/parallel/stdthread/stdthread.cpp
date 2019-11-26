// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file implements stdthread.hpp.

#include "dca/config/threading.hpp"

namespace dca {
namespace parallel {

constexpr char stdthread::parallel_type_str_[];

std::ostream& operator<<(std::ostream& o, const stdthread& c) {
  o << '\n'
    << "threading type:" << c.parallel_type_str_ << '\n'
    << "number of std::threads:" << ThreadPool::get_instance().size();
  return o;
}

}  // parallel
}  // dca
