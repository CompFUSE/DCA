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
// This class implements a fermionic Hilbert space.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_HILBERT_SPACES_HILBERT_SPACE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_HILBERT_SPACES_HILBERT_SPACE_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/basis_states/psi_operators.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/basis_states/psi_state.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/hilbert_spaces/hilbert_space_phi_representation.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

template <typename parameter_type, typename ed_options>  // N: size of bitset sequence
class Hilbert_space {
public:
  typedef typename ed_options::int_type int_type;

  typedef typename ed_options::scalar_type scalar_type;
  typedef typename ed_options::complex_type complex_type;

  typedef typename ed_options::phi_type phi_type;

public:
  Hilbert_space(const std::vector<std::string>& names_, const std::vector<int>& eigenvalues_);

  Hilbert_space(const Hilbert_space<parameter_type, ed_options>& HS);

  void insert(const psi_state<parameter_type, ed_options>& psi);

  void remove(int i) {
    psi_states.erase(psi_states.begin() + i);
  }
  bool remove(const psi_state<parameter_type, ed_options>& psi);

  bool mark_state(const psi_state<parameter_type, ed_options>& psi);

  typename std::vector<psi_state<parameter_type, ed_options>>::iterator binary_find(
      const psi_state<parameter_type, ed_options>& psi);

  void print(bool full = false);

  void sort_wrt_magnetization();
  void sort();

  int size() const {
    return psi_states.size();
  }
  int rep_size() const {
    return rep.size();
  }

  bool empty() {
    return psi_states.empty();
  }

  psi_state<parameter_type, ed_options>& get_element(int state) {
    return psi_states[state];
  }

  int get_N_occ() {
    return psi_states[0].occupation_number;
  }
  std::vector<std::string> get_name() {
    return names;
  }
  std::vector<int> get_eigenvalues() {
    return eigenvalues;
  }

  // assert or boolean to check wether occupation and magnetization is defined
  int get_occupation() {
    return eigenvalues[0];
  }
  int get_magnetization() {
    return eigenvalues[1];
  }

  Hilbert_space_phi_representation<parameter_type, ed_options>& get_rep() {
    return rep;
  }

  typename std::vector<psi_state<parameter_type, ed_options>>::iterator psi_states_begin() {
    return psi_states.begin();
  }

  typename std::vector<psi_state<parameter_type, ed_options>>::iterator psi_states_end() {
    return psi_states.end();
  }

  void initialize_rep() {
    rep.initialize(*this);
  }

private:
  std::vector<std::string> names;
  std::vector<int> eigenvalues;

  std::vector<psi_state<parameter_type, ed_options>> psi_states;

  Hilbert_space_phi_representation<parameter_type, ed_options> rep;
};

template <typename parameter_type, typename ed_options>
Hilbert_space<parameter_type, ed_options>::Hilbert_space(const std::vector<std::string>& names_,
                                                         const std::vector<int>& eigenvalues_)
    : names(names_), eigenvalues(eigenvalues_) {}

template <typename parameter_type, typename ed_options>
Hilbert_space<parameter_type, ed_options>::Hilbert_space(
    const Hilbert_space<parameter_type, ed_options>& HS)
    : names(HS.names), eigenvalues(HS.eigenvalues), psi_states(HS.psi_states) {}

template <typename parameter_type, typename ed_options>
void Hilbert_space<parameter_type, ed_options>::insert(const psi_state<parameter_type, ed_options>& psi) {
  psi_states.push_back(psi);
}

template <typename parameter_type, typename ed_options>
bool Hilbert_space<parameter_type, ed_options>::remove(const psi_state<parameter_type, ed_options>& psi) {
  typename std::vector<psi_state<parameter_type, ed_options>>::iterator it =
      std::find(psi_states.begin(), psi_states.end(), psi);

  if (it != psi_states.end()) {
    psi_states.erase(it);
    return true;
  }

  return false;
}

template <typename parameter_type, typename ed_options>
void Hilbert_space<parameter_type, ed_options>::print(bool full) {
  std::cout << size() << "\t" << rep.size();
  for (int_type i = 0; i < eigenvalues.size(); ++i) {
    std::cout << "\t" << eigenvalues[i];
  }
  std::cout << "\n";
  if (full) {
    for (int_type i = 0; i < psi_states.size(); ++i) {
      psi_states[i].print();
      std::cout << "\n";
    }
  }
}

template <typename parameter_type, typename ed_options>
void Hilbert_space<parameter_type, ed_options>::sort_wrt_magnetization() {
  std::sort(psi_states.begin(), psi_states.end(), less_magnetization<parameter_type, ed_options>);
}

template <typename parameter_type, typename ed_options>
void Hilbert_space<parameter_type, ed_options>::sort() {
  std::sort(psi_states.begin(), psi_states.end());
}

template <typename parameter_type, typename ed_options>
bool Hilbert_space<parameter_type, ed_options>::mark_state(
    const psi_state<parameter_type, ed_options>& psi) {
  typename std::vector<psi_state<parameter_type, ed_options>>::iterator it = binary_find(psi);

  if (it != psi_states.end()) {
    it->mark();
    return true;
  }
  else
    return false;
}

// Make sure that the subspace is sorted
template <typename parameter_type, typename ed_options>
typename std::vector<psi_state<parameter_type, ed_options>>::iterator Hilbert_space<
    parameter_type, ed_options>::binary_find(const psi_state<parameter_type, ed_options>& psi) {
  typename std::vector<psi_state<parameter_type, ed_options>>::iterator it =
      std::lower_bound(psi_states.begin(), psi_states.end(), psi);

  while (it != psi_states.end() &&
         !(ed::operator< <parameter_type, ed_options>(psi, *it))) {
    if (std::abs(abs(scalar_product<parameter_type, ed_options>(*it, psi)) - scalar_type(1.)) <
        ed_options::get_epsilon())  // 1.e-10)
      return it;
    ++it;
  }

  return psi_states.end();
}

}  // ed
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_HILBERT_SPACES_HILBERT_SPACE_HPP
