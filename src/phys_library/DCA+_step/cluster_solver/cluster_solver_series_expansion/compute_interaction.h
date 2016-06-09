// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// This class implements the interaction matrix.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_COMPUTE_INTERACTION_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_COMPUTE_INTERACTION_H

#include "comp_library/function_library/include_function_library.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"

namespace DCA {
namespace SERIES_EXPANSION {

class compute_interaction {
public:
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index

  using function_type = FUNC_LIB::function<double, dmn_2<nu, nu>>;

public:
  compute_interaction() {}

  template <class r_dmn_t>
  void execute(FUNC_LIB::function<double, dmn_3<nu, nu, r_dmn_t>>& H_interation) {
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
  FUNC_LIB::function<double, dmn_2<nu, nu>> U;
};
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_COMPUTE_INTERACTION_H
