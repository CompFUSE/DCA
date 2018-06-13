// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Bart Ydens
//
// This class organizes a vertex, i.e. it stores the time points of a \f$c\f$ and \f$c^{\dagger}\f$.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_STRUCTURES_HYBRIDIZATION_VERTEX_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_STRUCTURES_HYBRIDIZATION_VERTEX_HPP

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

class Hybridization_vertex {
public:
  typedef Hybridization_vertex this_type;

  Hybridization_vertex() : t_start_val(0), t_end_val(0) {}
  Hybridization_vertex(double t_start, double t_end) : t_start_val(t_start), t_end_val(t_end) {}

  double t_start() const {
    return t_start_val;
  }
  double t_end() const {
    return t_end_val;
  }

  void set_t_start(double t_start) {
    t_start_val = t_start;
  }
  void set_t_end(double t_end) {
    t_end_val = t_end;
  }

  this_type& operator=(this_type& other_vertex) {
    t_start_val = other_vertex.t_start();
    t_end_val = other_vertex.t_end();
    return *this;
  }
  this_type& operator=(const this_type& other_vertex) {
    t_start_val = other_vertex.t_start();
    t_end_val = other_vertex.t_end();
    return *this;
  }

private:
  double t_start_val, t_end_val;
};

inline bool operator<(const Hybridization_vertex& t1, const Hybridization_vertex& t2) {
  return t1.t_start() < t2.t_start();
}

inline bool operator<(const Hybridization_vertex& t1, const double t2) {
  return t1.t_start() < t2;
}

inline bool operator>(Hybridization_vertex t1, Hybridization_vertex t2) {
  return t1.t_start() > t2.t_start();
}

inline bool operator==(Hybridization_vertex t1, Hybridization_vertex t2) {
  return t1.t_start() == t2.t_start();
}

}  // cthyb
}  // phys
}  // solver
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_STRUCTURES_HYBRIDIZATION_VERTEX_HPP
