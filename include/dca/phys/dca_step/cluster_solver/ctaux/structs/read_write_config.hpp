// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the method to read and write the CT-AUX configuration.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_READ_WRITE_CONFIG_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_READ_WRITE_CONFIG_HPP

#include "dca/io/buffer.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/ct_aux_hs_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_pair.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <class Parameters>
io::Buffer& operator<<(io::Buffer& buff, const vertex_pair<Parameters>& v) {
  return buff << v.bands << v.e_spins << v.spin_orbitals << v.r_sites << v.HS_spin << v.tau;
}

template <class Parameters>
io::Buffer& operator>>(io::Buffer& buff, vertex_pair<Parameters>& v) {
  v.creatable = false;
  v.annihilatable = true;
  v.shuffled = true;

  return buff >> v.bands >> v.e_spins >> v.spin_orbitals >> v.r_sites >> v.HS_spin >> v.tau;
}

template <class Parameters>
io::Buffer& operator<<(io::Buffer& buff, const CT_AUX_HS_configuration<Parameters>& config) {
  buff << config.configuration.size();
  for (const auto& vertex : config.configuration)
    buff << vertex;
  return buff;
}

template <class Parameters>
io::Buffer& operator>>(io::Buffer& buff, CT_AUX_HS_configuration<Parameters>& config) {
  config.reset();

  std::size_t n;
  buff >> n;

  for (int i = 0; i < n; ++i) {
    vertex_pair<Parameters> vertex(config.parameters, config.rng, config.configuration.size(),
                                   config.configuration_e_DN.size(),
                                   config.configuration_e_UP.size(), config.next_vertex_id_++);

    buff >> vertex;

    ++config.current_Nb_of_annihilatable_spins;
    config.update_configuration_e_spin(vertex);
    config.configuration.push_back(vertex);
  }
  return buff;
}

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_READ_WRITE_CONFIG_HPP
