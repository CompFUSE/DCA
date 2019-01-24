// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements vertex_singleton.hpp

#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"

#include <cassert>
#include <iostream>

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

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

bool vertex_singleton::operator==(const dca::phys::solver::ctaux::vertex_singleton& rhs) const {
  return band == rhs.get_band() && e_spin == rhs.get_e_spin() &&
         spin_orbital == rhs.get_spin_orbital()

         && paired_spin_orbital == rhs.get_paired_spin_orbital() && r_site == rhs.get_r_site() &&
         delta_r == rhs.get_delta_r() && tau == rhs.get_tau()

         && HS_spin == rhs.get_HS_spin() && HS_field == rhs.get_HS_field() &&
         configuration_index == rhs.get_configuration_index();

  //  if (!result) {
  //    std::cout << band << "\t" << rhs.get_band() << std::endl;
  //    std::cout << e_spin << "\t" << rhs.get_e_spin() << std::endl;
  //    std::cout << spin_orbital << "\t" << rhs.get_spin_orbital() << std::endl;
  //
  //    std::cout << HS_field << "\t" << rhs.get_HS_field() << std::endl;
  //  }
  //  return result;
}

}  // ctaux
}  // solver
}  // phys
}  // dca
