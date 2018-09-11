// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements no_concurrency.hpp.

#include "dca/parallel/no_concurrency/no_concurrency.hpp"

#include <iostream>

namespace dca {
namespace parallel {

void NoConcurrency::abort() const {
  std::cout << "\nAborting process.\n";
  std::terminate();
}

constexpr char NoConcurrency::parallel_type_str_[];

std::ostream& operator<<(std::ostream& o, const NoConcurrency& c) {
  o << '\n'
    << "concurrency type:" << c.parallel_type_str_ << '\n'
    << "number of processors:" << c.number_of_processors() << '\n'
    << "grouping first:" << c.first() << '\n'
    << "grouping last::" << c.last();
  return o;
}

}  // parallel
}  // dca
