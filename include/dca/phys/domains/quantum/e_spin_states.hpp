// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
