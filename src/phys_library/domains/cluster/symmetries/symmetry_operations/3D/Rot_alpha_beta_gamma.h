// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Elements of the cyclic group C_n(m).

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_ROT_ALPHA_BETA_GAMMA_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_ROT_ALPHA_BETA_GAMMA_H

#include <cmath>
#include "phys_library/domains/cluster/symmetries/symmetry_operations/group_action.h"

template <int alpha_num, int alpha_den, int beta_num, int beta_den, int gamma_num, int gamma_den>
class rot_a_b_c : public group_action<3> {
public:
  typedef group_action<3> base_type;
  typedef rot_a_b_c<alpha_num, alpha_den, beta_num, beta_den, gamma_num, gamma_den> this_type;

  rot_a_b_c(){};

  const static double* matrix() {
    double a = 2 * M_PI * double(alpha_num) / double(alpha_den);  // polar coordinate \theta
    double b = 2 * M_PI * double(beta_num) / double(beta_den);    // polar coordinate \phi
    double c = 2 * M_PI * double(gamma_num) / double(gamma_den);  // actual rotation angle

    // matrix in column-major format...
    static double matrix[3 * 3] = {
        cos(a) * (cos(a) * cos(c) - cos(b) * sin(a) * sin(c)) -
            sin(a) * (-(sin(a) * pow(sin(b), 2)) +
                      cos(b) * (-(cos(b) * cos(c) * sin(a)) - cos(a) * sin(c))),
        cos(a) * (cos(c) * sin(a) + cos(a) * cos(b) * sin(c)) -
            sin(a) * (cos(a) * pow(sin(b), 2) + cos(b) * (cos(a) * cos(b) * cos(c) - sin(a) * sin(c))),
        -(sin(a) * (-(cos(b) * sin(b)) + cos(b) * cos(c) * sin(b))) + cos(a) * sin(b) * sin(c),

        sin(a) * (cos(a) * cos(c) - cos(b) * sin(a) * sin(c)) +
            cos(a) * (-(sin(a) * pow(sin(b), 2)) +
                      cos(b) * (-(cos(b) * cos(c) * sin(a)) - cos(a) * sin(c))),
        sin(a) * (cos(c) * sin(a) + cos(a) * cos(b) * sin(c)) +
            cos(a) * (cos(a) * pow(sin(b), 2) + cos(b) * (cos(a) * cos(b) * cos(c) - sin(a) * sin(c))),
        cos(a) * (-(cos(b) * sin(b)) + cos(b) * cos(c) * sin(b)) + sin(a) * sin(b) * sin(c),

        cos(b) * sin(a) * sin(b) + sin(b) * (-(cos(b) * cos(c) * sin(a)) - cos(a) * sin(c)),
        -(cos(a) * cos(b) * sin(b)) + sin(b) * (cos(a) * cos(b) * cos(c) - sin(a) * sin(c)),
        pow(cos(b), 2) + cos(c) * pow(sin(b), 2)};
    return matrix;
  }

  /*
  // MATHEMATICA CODE

  // Z[\[Alpha]_] := ({
  //    {Cos[\[Alpha]], -Sin[\[Alpha]], 0},
  //    {Sin[\[Alpha]], Cos[\[Alpha]], 0},
  //    {0, 0, 1}
  //   })
  // X[\[Alpha]_] := ({
  //    {1, 0, 0},
  //    {0, Cos[\[Alpha]], -Sin[\[Alpha]]},
  //    {0, Sin[\[Alpha]], Cos[\[Alpha]]}
  //   })
  //
  // Transpose[
  //   Z[\[Alpha]].X[\[Beta]].Z[\[Gamma]].X[-\[Beta]].Z[-\[Alpha]]]
  // matrix = CForm[Transpose[
  //    X[\[Alpha]].Z[\[Beta]].X[\[Gamma]].Z[-\[Beta]].X[-\[Alpha]]] /. {\
  // \[Alpha] -> a, \[Beta] -> b, \[Gamma] -> g}]
  */
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_ROT_ALPHA_BETA_GAMMA_H
