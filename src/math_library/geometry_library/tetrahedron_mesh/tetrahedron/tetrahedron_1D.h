// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Full template specialization for a 1D tetrahedron.

#ifndef MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_TETRAHEDRON_1D_H
#define MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_TETRAHEDRON_1D_H

#include "math_library/geometry_library/tetrahedron_mesh/tetrahedron/tetrahedron_template.h"

#include <cassert>
#include <cstdlib>
#include <vector>

#include "comp_library/function_plotting/include_plotting.h"

namespace math_algorithms {

template <>
struct tetrahedron<1> {
public:
  tetrahedron();
  ~tetrahedron();

  std::vector<double> compute_cm();

  void translate(std::vector<double> q);

  void plot(Gnuplot& plot_obj);
  void plot_q_vecs(Gnuplot& plot_obj);

  template <typename scalartype>
  void update_gaussian_domain(int& size, std::vector<scalartype>& weights,
                              std::vector<std::vector<scalartype>>& k_vecs);

  template <typename scalartype>
  void update_tetrahedron_domain(int& size, std::vector<scalartype>& weights,
                                 std::vector<std::vector<scalartype>>& k_vecs);

  template <typename mesh_t>
  void do_recursion(std::vector<tetrahedron<1>>& tetrahedra, mesh_t& mesh);

public:
  int index[2];

  double volume;

  std::vector<double> vec_0;
  std::vector<double> vec_1;

  std::vector<double> normal;

  int N_q;

  double* q_w;
  double* q_vecs;
};

tetrahedron<1>::tetrahedron() {
  N_q = -1;

  q_w = NULL;
  q_vecs = NULL;
}

tetrahedron<1>::~tetrahedron() {
  if (q_w != NULL)
    delete[] q_w;

  if (q_vecs != NULL)
    delete[] q_vecs;
}

std::vector<double> tetrahedron<1>::compute_cm() {
  std::vector<double> cm(1, 0);

  cm[0] = (vec_0[0] + vec_1[0]) / 2.;

  return cm;
}

void tetrahedron<1>::translate(std::vector<double> k_point) {
  assert(k_point.size() == 1);

  for (int j = 0; j < 1; j++) {
    vec_0[j] += k_point[j];
    vec_1[j] += k_point[j];
  }

  for (int i = 0; i < N_q; i++)
    for (int j = 0; j < 1; j++)
      q_vecs[j + i * 1] += k_point[j];
}

void tetrahedron<1>::plot(Gnuplot& plot_obj) {
  std::vector<double> x(2);
  std::vector<double> y(2);

  x[0] = vec_0[0];
  y[0] = 0;
  x[1] = vec_1[0];
  y[1] = 0;

  plot_obj.plot_xy(x, y);
}

void tetrahedron<1>::plot_q_vecs(Gnuplot& plot_obj) {
  std::vector<double> x(N_q);
  std::vector<double> y(N_q);

  for (int i = 0; i < N_q; i++) {
    x[i] = q_vecs[i];
    y[i] = 0;
  }

  plot_obj.plot_xy(x, y);
}

template <typename mesh_t>
void tetrahedron<1>::do_recursion(std::vector<tetrahedron<1>>& tetrahedra, mesh_t& mesh) {
  std::vector<double> k0 = mesh[index[0]];
  std::vector<double> k1 = mesh[index[1]];

  static std::vector<double> km(1, 0);
  km[0] = (k0[0] + k1[0]) / 2.;

  int ind = mesh.size();
  mesh.push_back(km);

  static tetrahedron<1> tet_l;
  {
    tet_l.index[0] = index[0];
    tet_l.index[1] = ind;

    tet_l.vec_0 = k0;
    tet_l.vec_1 = km;

    tet_l.normal = this->normal;
    tet_l.volume = std::abs(km[0] - k0[0]);
  }

  static tetrahedron<1> tet_r;
  {
    tet_r.index[0] = ind;
    tet_r.index[1] = index[1];

    tet_r.vec_0 = km;
    tet_r.vec_1 = k1;

    tet_r.volume = std::abs(k1[0] - km[0]);
    tet_r.normal = normal;
  }

  tetrahedra.push_back(tet_l);
  tetrahedra.push_back(tet_r);
}
}

#endif  // MATH_LIBRARY_GEOMETRY_LIBRARY_TETRAHEDRON_MESH_TETRAHEDRON_TETRAHEDRON_1D_H
