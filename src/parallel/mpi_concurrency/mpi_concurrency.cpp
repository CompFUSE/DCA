// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file implements mpi_concurrency.hpp.

#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"

namespace dca {
namespace parallel {

constexpr char MPIConcurrency::parallel_type_str_[];

std::ostream& operator<<(std::ostream& o, const MPIConcurrency& c) {
  o << '\n'
    << "concurrency type:" << c.parallel_type_str_ << '\n'
    << "number of processors:" << c.number_of_processors() << '\n'
    << "grouping first:" << c.first() << '\n'
    << "grouping last::" << c.last();
  return o;
}

}  // parallel
}  // dca
