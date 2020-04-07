// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author:Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This class organizes the vertex configuration for ct-int.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_UTIL_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_UTIL_HPP

#include <stdexcept>
#include <vector>

#include "dca/util/random_access_map.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

using VertexTypeList = dca::util::RandomAccessMap<std::size_t, std::size_t, 16>;

std::size_t getRandomElement(const std::vector<const VertexTypeList*>& v_ptrs,
                          double rand) noexcept;

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_UTIL_HPP
