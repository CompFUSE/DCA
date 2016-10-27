// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides creation and annihilation operators for phi states.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_PHI_FERMIONIC_OPERATORS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_PHI_FERMIONIC_OPERATORS_HPP

#include <cassert>

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

template <typename parameter_type, typename ed_options>  // N: size of bitset sequence
class PhiFermionicOperators {
public:
  typedef typename ed_options::int_type int_type;
  typedef typename ed_options::phi_type phi_type;

  static bool create_at(int_type l, phi_type& phi, int& sign);
  static bool annihilate_at(int_type l, phi_type& phi, int& sign);
};

template <typename parameter_type, typename ed_options>
bool PhiFermionicOperators<parameter_type, ed_options>::create_at(int_type l, phi_type& phi,
                                                                  int& sign) {
  int tmp_sgn = sign;

  if (phi.test(l)) {
    return false;
  }
  else {
    phi.set(l);

    if (l != 0) {
      int_type tmp = (1 << l) - 1;

      phi_type mask(tmp);

      bool change_sign = ((phi & mask).count()) & 1;

      if (change_sign)
        sign *= -1;
    }

    for (int i = 0; i < l; ++i) {
      if (phi[i])
        tmp_sgn *= -1;
    }

    assert(tmp_sgn == sign);

    return true;
  }
}

template <typename parameter_type, typename ed_options>
bool PhiFermionicOperators<parameter_type, ed_options>::annihilate_at(int_type l, phi_type& phi,
                                                                      int& sign) {
  int tmp_sgn = sign;

  if (phi.test(l)) {
    phi.reset(l);

    if (l != 0) {
      int_type tmp = (1 << l) - 1;

      phi_type mask(tmp);

      bool change_sign = ((phi & mask).count()) & 1;

      if (change_sign)
        sign *= -1;
    }

    for (int i = 0; i < l; ++i) {
      if (phi[i])
        tmp_sgn *= -1;
    }

    assert(tmp_sgn == sign);

    return true;
  }
  else {
    return false;
  }
}

}  // ed
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_PHI_FERMIONIC_OPERATORS_HPP
