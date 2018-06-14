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

// Modifies the state 'phi' by creating an electron in the l-th spin-orbital.
// Changes the overall sign of the state, if required by the fermionic algebra.
template <typename parameter_type, typename ed_options>
bool PhiFermionicOperators<parameter_type, ed_options>::create_at(int_type l, phi_type& phi,
                                                                  int& sign) {
  // Spin-orbital already occupied.
  if (phi.test(l)) {
    return false;
  }

  else {
    // Create electron in l-th spin-orbital.
    phi.set(l);

// Compute the new sign the easy (slow) way.
#ifndef NDEBUG
    int new_sgn = sign;
    for (int i = 0; i < l; ++i) {
      if (phi[i])
        new_sgn *= -1;
    }
#endif  // NDEBUG

    // Determine whether there is a sign change in the elegant (fast) way.
    if (l != 0) {
      int_type tmp = (1 << l) - 1;
      phi_type mask(tmp);  // = 00...0011...11,
                           // where the first 0 (from the right) is at l-th position.

      // Determines the parity of the number of 1's in phi before the l-th position:
      // odd --> change_sign = true
      // even --> change_sign = false.
      bool change_sign = ((phi & mask).count()) & 1;

      if (change_sign)
        sign *= -1;
    }

    assert(new_sgn == sign);

    return true;
  }
}

// Modifies the state 'phi' by annihilating the electron in the l-th spin-orbital.
// Changes the overall sign of the state, if required by the fermionic algebra.
template <typename parameter_type, typename ed_options>
bool PhiFermionicOperators<parameter_type, ed_options>::annihilate_at(int_type l, phi_type& phi,
                                                                      int& sign) {
  if (phi.test(l)) {
    // Annihilate the electron in the l-th spin-oribital.
    phi.reset(l);

// Compute the new sign the easy (slow) way.
#ifndef NDEBUG
    int new_sgn = sign;
    for (int i = 0; i < l; ++i) {
      if (phi[i])
        new_sgn *= -1;
    }
#endif  // NDEBUG

    // Determine whether there is a sign change in the elegant (fast) way.
    if (l != 0) {
      int_type tmp = (1 << l) - 1;
      phi_type mask(tmp);  // = 00...0011...11,
                           // where the first 0 (from the right) is at l-th position.

      // Determines the parity of the number of 1's in phi before the l-th position:
      // odd --> change_sign = true
      // even --> change_sign = false.
      bool change_sign = ((phi & mask).count()) & 1;

      if (change_sign)
        sign *= -1;
    }

    assert(new_sgn == sign);

    return true;
  }

  // Spin-orbital not occupied.
  else {
    return false;
  }
}

}  // ed
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_PHI_FERMIONIC_OPERATORS_HPP
