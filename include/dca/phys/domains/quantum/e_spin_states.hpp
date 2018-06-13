// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the two electron spin states, up and down.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_E_SPIN_STATES_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_E_SPIN_STATES_HPP

namespace dca {
namespace phys {
// dca::phys::

enum e_spin_states { e_DN = -1, e_UP = 1 };
using e_spin_states_type = e_spin_states;

}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_E_SPIN_STATES_HPP
