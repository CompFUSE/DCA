// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//         Bart Ydens
//
// This class organizes a vertex, stores the time points of a \f$c\f$ and \f$c^{\dagger}\f$.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_STRUCTURES_SS_HYBRIDIZATION_VERTEX_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_STRUCTURES_SS_HYBRIDIZATION_VERTEX_H

namespace DCA {

class Hybridization_vertex {
public:
  typedef Hybridization_vertex this_type;

public:
  Hybridization_vertex();
  Hybridization_vertex(double t_start, double t_end);

  double t_start() const;
  double t_end() const;

  void set_t_start(double t_start);
  void set_t_end(double t_end);

  this_type& operator=(this_type& other_vertex);
  this_type& operator=(const this_type& other_vertex);

private:
  double t_start_val, t_end_val;
};

Hybridization_vertex::Hybridization_vertex() : t_start_val(0), t_end_val(0) {}

Hybridization_vertex::Hybridization_vertex(double t_st, double t_e)
    : t_start_val(t_st), t_end_val(t_e) {}

double Hybridization_vertex::t_start() const {
  return t_start_val;
}

double Hybridization_vertex::t_end() const {
  return t_end_val;
}

void Hybridization_vertex::set_t_start(double t_st) {
  t_start_val = t_st;
}

void Hybridization_vertex::set_t_end(double t_e) {
  t_end_val = t_e;
}

bool operator<(const Hybridization_vertex& t1, const Hybridization_vertex& t2) {
  return t1.t_start() < t2.t_start();
}

bool operator<(const Hybridization_vertex& t1, const double t2) {
  return t1.t_start() < t2;
}

bool operator>(Hybridization_vertex t1, Hybridization_vertex t2) {
  return t1.t_start() > t2.t_start();
}

bool operator==(Hybridization_vertex t1, Hybridization_vertex t2) {
  return t1.t_start() == t2.t_start();
}

Hybridization_vertex& Hybridization_vertex::operator=(Hybridization_vertex& other_vertex) {
  t_start_val = other_vertex.t_start();
  t_end_val = other_vertex.t_end();

  return *this;
}

Hybridization_vertex& Hybridization_vertex::operator=(const Hybridization_vertex& other_vertex) {
  t_start_val = other_vertex.t_start();
  t_end_val = other_vertex.t_end();

  return *this;
}

}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SS_HYBRIDIZATION_SS_HYBRIDIZATION_STRUCTURES_SS_HYBRIDIZATION_VERTEX_H
