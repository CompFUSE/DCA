// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_FUNCTION_PROXY_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_FUNCTION_PROXY_HPP

#include "dca/function/function.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::pyhs::solver::ctint::

template<typename ScalarType, class Domain>
class FunctionProxy {
  public:
  template<class OtherDmn>
  FunctionProxy(const func::function<ScalarType, OtherDmn>& f):
  fnc_values_(f.values()),
  dmn_(){
    if(dmn_.get_size() != f.get_domain().get_size())
      throw(std::logic_error("Domain sizes are different.\n"));
  }

  template <typename... Ts>
  ScalarType operator()(Ts... indices) const{
    return fnc_values_[dmn_(indices...)];
  }

private:
  const ScalarType* fnc_values_;
  Domain dmn_;
};

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_FUNCTION_PROXY_HPP
