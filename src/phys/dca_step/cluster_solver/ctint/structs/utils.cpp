// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author:Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// Helper methods implementation fot the ct-int configuration.

#include "dca/phys/dca_step/cluster_solver/ctint/structs/utils.hpp"

#include <cassert>
#include <vector>

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

unsigned getRandomElement(const std::vector<const std::vector<std::size_t>*>& v_ptrs,
                          const double rand) noexcept {
  assert(rand >= 0 && rand <= 1);

  unsigned size = 0;
  for (int i = v_ptrs.size() - 1; i >= 0; --i)
    size += v_ptrs[i]->size();

  if (size == 0)
    return -1;

  // TODO: use binary search or other efficient scheme.
  unsigned idx = size * rand;
  for (auto v_ptr : v_ptrs) {
    if (idx < v_ptr->size())
      return (*v_ptr)[idx];
    idx -= v_ptr->size();
  }

  return -1;
}

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
