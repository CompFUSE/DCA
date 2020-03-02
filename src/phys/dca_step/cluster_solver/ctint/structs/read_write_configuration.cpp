// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author:Giovanni Balduzzi (gbalduzz@phys.ethz.ch)
//
// This file defines the I/O of the CT-INT configuration.

#include "dca/phys/dca_step/cluster_solver/ctint/structs/solver_configuration.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

io::Buffer& operator<<(io::Buffer& buff, const Vertex& v) {
  return buff << v.aux_spin << v.interaction_id << v.tau;
}

io::Buffer& operator>>(io::Buffer& buff, Vertex& v) {
  return buff >> v.aux_spin >> v.interaction_id >> v.tau;
}

io::Buffer& operator<<(io::Buffer& buff, const SolverConfiguration& config) {
  buff << config.vertices_.size();

  for (const auto& v : config.vertices_)
    buff << v;

  return buff;
}

io::Buffer& operator>>(io::Buffer& buff, SolverConfiguration& config) {
  std::size_t n;
  buff >> n;

  for (int i = 0; i < n; ++i) {
    Vertex v;
    buff >> v;
    v.tag = config.current_tag_++;
    config.push_back(v);
    config.commitInsertion(config.size() - 1);
  }

  return buff;
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
