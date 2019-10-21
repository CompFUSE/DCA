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

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

template<class T>
T getRandomElement(const std::vector<const std::vector<T>*>& v_ptrs, const double rand){
    assert(rand >= 0 && rand <= 1);

    unsigned size = 0;
    for(int i = v_ptrs.size() - 1; i >= 0; --i)
        size += v_ptrs[i]->size();

    // TODO: use binary search or other efficient scheme.
    unsigned idx = size * rand;
    for(auto v_ptr : v_ptrs){
        if (idx < v_ptr->size())
            return (*v_ptr)[idx];
        idx -= v_ptr->size();
    }

    throw(std::logic_error("Something went wrong."));
}


}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_STRUCTS_UTIL_HPP
