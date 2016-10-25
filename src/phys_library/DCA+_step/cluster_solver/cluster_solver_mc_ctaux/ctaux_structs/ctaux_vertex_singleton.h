// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class represents a vertex singleton.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_STRUCTS_CTAUX_VERTEX_SINGLETON_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_STRUCTS_CTAUX_VERTEX_SINGLETON_H

#include <cassert>
#include <iostream>

#include "dca/phys/domains/quantum/e_spin_states.hpp"

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_domains/HS_field_sign_domain.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_domains/HS_spin_domain.h"

using namespace dca::phys;

namespace DCA {
namespace QMCI {
// DCA::QMCI::

class vertex_singleton {
public:
  typedef vertex_singleton this_type;

public:
  vertex_singleton();

  //       vertex_singleton(int    band_in,
  // 		       int    r_site_in,
  // 		       double tau_in);

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

  int get_band() const;
  e_spin_states_type get_e_spin() const;
  int get_spin_orbital() const;
  int get_paired_spin_orbital() const;
  int get_r_site() const;
  int get_delta_r() const;
  double get_tau() const;
  HS_spin_states_type get_HS_spin() const;
  HS_field_sign_type get_HS_field() const;
  int get_configuration_index() const;

  HS_spin_states_type& get_HS_spin();
  int& get_configuration_index();

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

vertex_singleton::vertex_singleton() {}

//     vertex_singleton::vertex_singleton(int    band_in,
// 				       int    r_site_in,
// 				       double tau_in):
//       band(band_in),
//       r_site(r_site_in),
//       tau(tau_in)
//     {}

vertex_singleton::vertex_singleton(int band_in, e_spin_states_type e_spin_in, int spin_orbital_in,

                                   int paired_spin_orbital_in, int r_site_in, int delta_r_in,
                                   double tau_in,

                                   HS_spin_states_type HS_spin_in, HS_field_sign_type HS_field_in,
                                   int configuration_index_in)
    : band(band_in),
      e_spin(e_spin_in),
      spin_orbital(spin_orbital_in),

      paired_spin_orbital(paired_spin_orbital_in),
      r_site(r_site_in),
      delta_r(delta_r_in),
      tau(tau_in),

      HS_spin(HS_spin_in),
      HS_field(HS_field_in),
      configuration_index(configuration_index_in) {}

vertex_singleton::vertex_singleton(const vertex_singleton& other_vertex_couple)
    : band(other_vertex_couple.get_band()),
      e_spin(other_vertex_couple.get_e_spin()),
      spin_orbital(other_vertex_couple.get_spin_orbital()),

      paired_spin_orbital(other_vertex_couple.get_paired_spin_orbital()),
      r_site(other_vertex_couple.get_r_site()),
      delta_r(other_vertex_couple.get_delta_r()),
      tau(other_vertex_couple.get_tau()),

      HS_spin(other_vertex_couple.get_HS_spin()),
      HS_field(other_vertex_couple.get_HS_field()),
      configuration_index(other_vertex_couple.get_configuration_index()) {}

vertex_singleton& vertex_singleton::operator=(vertex_singleton& other_vertex_couple) {
  band = other_vertex_couple.get_band();
  e_spin = other_vertex_couple.get_e_spin();
  spin_orbital = other_vertex_couple.get_spin_orbital();

  paired_spin_orbital = other_vertex_couple.get_paired_spin_orbital();
  r_site = other_vertex_couple.get_r_site();
  delta_r = other_vertex_couple.get_delta_r();
  tau = other_vertex_couple.get_tau();

  HS_spin = other_vertex_couple.get_HS_spin();
  HS_field = other_vertex_couple.get_HS_field();
  configuration_index = other_vertex_couple.get_configuration_index();

  return *this;
}

vertex_singleton& vertex_singleton::operator=(const vertex_singleton& other_vertex_couple) {
  band = other_vertex_couple.get_band();
  e_spin = other_vertex_couple.get_e_spin();
  spin_orbital = other_vertex_couple.get_spin_orbital();

  paired_spin_orbital = other_vertex_couple.get_paired_spin_orbital();
  r_site = other_vertex_couple.get_r_site();
  delta_r = other_vertex_couple.get_delta_r();
  tau = other_vertex_couple.get_tau();

  HS_spin = other_vertex_couple.get_HS_spin();
  HS_field = other_vertex_couple.get_HS_field();
  configuration_index = other_vertex_couple.get_configuration_index();

  return *this;
}

bool vertex_singleton::equals(vertex_singleton other_vertex_couple) {
  bool result;

  if (band == other_vertex_couple.get_band() && e_spin == other_vertex_couple.get_e_spin() &&
      spin_orbital == other_vertex_couple.get_spin_orbital()

      && paired_spin_orbital == other_vertex_couple.get_paired_spin_orbital() &&
      r_site == other_vertex_couple.get_r_site() && delta_r == other_vertex_couple.get_delta_r() &&
      tau == other_vertex_couple.get_tau()

      && HS_spin == other_vertex_couple.get_HS_spin() &&
      HS_field == other_vertex_couple.get_HS_field() &&
      configuration_index == other_vertex_couple.get_configuration_index())
    result = true;
  else
    result = false;

  if (!result) {
    std::cout << band << "\t" << other_vertex_couple.get_band() << std::endl;
    std::cout << e_spin << "\t" << other_vertex_couple.get_e_spin() << std::endl;
    std::cout << spin_orbital << "\t" << other_vertex_couple.get_spin_orbital() << std::endl;

    std::cout << HS_field << "\t" << other_vertex_couple.get_HS_field() << std::endl;
  }

  assert(result);

  return result;
}

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

int vertex_singleton::get_band() const {
  return band;
}

e_spin_states_type vertex_singleton::get_e_spin() const {
  return e_spin;
}

int vertex_singleton::get_spin_orbital() const {
  return spin_orbital;
}

int vertex_singleton::get_paired_spin_orbital() const {
  return paired_spin_orbital;
}

int vertex_singleton::get_r_site() const {
  return r_site;
}

int vertex_singleton::get_delta_r() const {
  return delta_r;
}

double vertex_singleton::get_tau() const {
  return tau;
}

HS_spin_states_type vertex_singleton::get_HS_spin() const {
  return HS_spin;
}

HS_spin_states_type& vertex_singleton::get_HS_spin() {
  return HS_spin;
}

HS_field_sign_type vertex_singleton::get_HS_field() const {
  return HS_field;
}

int vertex_singleton::get_configuration_index() const {
  return configuration_index;
}

int& vertex_singleton::get_configuration_index() {
  return configuration_index;
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_STRUCTS_CTAUX_VERTEX_SINGLETON_H
