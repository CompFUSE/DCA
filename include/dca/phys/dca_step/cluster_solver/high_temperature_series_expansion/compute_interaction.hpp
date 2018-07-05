// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class stores the interaction matrix.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_COMPUTE_INTERACTION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_COMPUTE_INTERACTION_HPP

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace htseries {
// dca::phys::solver::htseries::

class compute_interaction {
public:
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using function_type = func::function<double, func::dmn_variadic<nu, nu>>;

public:
  compute_interaction() {}

  template <class r_dmn_t>
  void execute(func::function<double, func::dmn_variadic<nu, nu, r_dmn_t>>& H_interation) {
    for (int nu0 = 0; nu0 < 2 * b::dmn_size(); ++nu0)
      for (int nu1 = 0; nu1 < 2 * b::dmn_size(); ++nu1)
        U(nu0, nu1) = H_interation(nu0, nu1, 0);

    // 	for(int s_0=0; s_0<2*b::dmn_size(); ++s_0){
    // 	  for(int s_1=0; s_1<2*b::dmn_size(); ++s_1)
    // 	    cout << U(s_0, s_1) << "\t";
    // 	  cout << "\n";
    // 	}
    // 	cout << "\n";
  }

  function_type& get_function() {
    return U;
  }

protected:
  func::function<double, func::dmn_variadic<nu, nu>> U;
};

}  // htseries
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_COMPUTE_INTERACTION_HPP
