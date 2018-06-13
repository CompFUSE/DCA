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
// This file implements get_bounds.hpp.

#include "dca/parallel/util/get_bounds.hpp"
#include <cassert>

namespace dca {
namespace parallel {
namespace util {
// dca::parallel::util::

std::pair<int, int> getBounds(int id, int num_threads, const std::pair<int, int>& current_bounds) {
  assert(id >= 0 && id < num_threads);

  const int range = current_bounds.second - current_bounds.first;
  assert(range >= 0);

  if (num_threads < range)
    return std::make_pair((id * range) / num_threads + current_bounds.first,
                          ((id + 1) * range) / num_threads + current_bounds.first);

  else if (id < range)
    return std::make_pair(id + current_bounds.first, id + 1 + current_bounds.first);

  else
    return std::make_pair(-1 + current_bounds.first, -1 + current_bounds.first);
}

}  // util
}  // parallel
}  // dca
