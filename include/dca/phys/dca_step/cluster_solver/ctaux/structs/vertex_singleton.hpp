// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements a vertex singleton.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_VERTEX_SINGLETON_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_VERTEX_SINGLETON_HPP

#include "dca/phys/dca_step/cluster_solver/ctaux/domains/hs_field_sign_domain.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/domains/hs_spin_domain.hpp"
#include "dca/phys/domains/quantum/e_spin_states.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

class vertex_singleton {
public:
  typedef vertex_singleton this_type;

public:
  vertex_singleton() {}

  vertex_singleton(int band_in, e_spin_states_type e_spin_in, int spin_orbital_in,
                   int paired_spin_orbital_in, int r_site_in, int delta_r_in, double tau_in,
                   HS_spin_states_type HS_spin_in, HS_field_sign_type HS_field_in,
                   int configuration_index_in);

  vertex_singleton(const this_type& other_vertex_couple);

  this_type& operator=(this_type& other_vertex_couple);

  this_type& operator=(const this_type& other_vertex_couple);

  bool equals(
      this_type other_vertex_couple);  // --> needed for consistency-check in HS-configuration !!

  template <class configuration_type>
  vertex_singleton& get_partner(configuration_type& configuration);

  int get_band() const {
    return band;
  }
  e_spin_states_type get_e_spin() const {
    return e_spin;
  }
  int get_spin_orbital() const {
    return spin_orbital;
  }
  int get_paired_spin_orbital() const {
    return paired_spin_orbital;
  }
  int get_r_site() const {
    return r_site;
  }
  int get_delta_r() const {
    return delta_r;
  }
  double get_tau() const {
    return tau;
  }
  HS_spin_states_type get_HS_spin() const {
    return HS_spin;
  }
  HS_spin_states_type& get_HS_spin() {
    return HS_spin;
  }
  HS_field_sign_type get_HS_field() const {
    return HS_field;
  }
  int get_configuration_index() const {
    return configuration_index;
  }
  int& get_configuration_index() {
    return configuration_index;
  }

  // Interface for accumulator.
  int get_left_band() const {
    return get_band();
  }
  int get_right_band() const {
    return get_band();
  }
  int get_left_site() const {
    return get_r_site();
  }
  int get_right_site() const {
    return get_r_site();
  }

private:
  int band;
  e_spin_states_type e_spin;
  int spin_orbital;
  int paired_spin_orbital;

  int r_site;
  int delta_r;
  double tau;

  HS_spin_states_type HS_spin;
  HS_field_sign_type HS_field;

  int configuration_index;
};

template <class configuration_type>
vertex_singleton& vertex_singleton::get_partner(configuration_type& configuration) {
  e_spin_states_type e_spin;
  int configuration_e_spin;

  if (HS_field == HS_FIELD_DN) {
    e_spin = configuration[configuration_index].get_e_spins().second;
    configuration_e_spin =
        configuration[configuration_index].get_configuration_e_spin_indices().second;
  }
  else {
    e_spin = configuration[configuration_index].get_e_spins().first;
    configuration_e_spin =
        configuration[configuration_index].get_configuration_e_spin_indices().first;
  }

  return configuration.get(e_spin)[configuration_e_spin];
}

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_VERTEX_SINGLETON_HPP
