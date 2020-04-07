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

std::size_t getRandomElement(const std::vector<const VertexTypeList*>& container_ptrs,
                             const double rand) noexcept {
  assert(rand >= 0 && rand <= 1);

  std::size_t size = 0;
  for (auto c_ptr : container_ptrs)
    size += c_ptr->size();

  if (size == 0)
    return -1;

  // TODO: use binary search or other efficient scheme.
  std::size_t idx = size * rand;
  for (auto c_ptr : container_ptrs) {
    if (idx < c_ptr->size())
      return (*c_ptr)[idx];
    idx -= c_ptr->size();
  }

  return -1;
}

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
