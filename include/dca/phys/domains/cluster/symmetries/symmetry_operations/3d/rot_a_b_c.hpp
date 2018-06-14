// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Elements of the cyclic group C_n(m).
// Rotation alpha beta gamma.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_ROT_A_B_C_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_ROT_A_B_C_HPP

#include <cmath>
#include "dca/phys/domains/cluster/symmetries/symmetry_operations/group_action.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

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
        std::cos(a) * (std::cos(a) * std::cos(c) - std::cos(b) * std::sin(a) * std::sin(c)) -
            std::sin(a) * (-(std::sin(a) * std::pow(std::sin(b), 2)) +
                           std::cos(b) * (-(std::cos(b) * std::cos(c) * std::sin(a)) -
                                          std::cos(a) * std::sin(c))),
        std::cos(a) * (std::cos(c) * std::sin(a) + std::cos(a) * std::cos(b) * std::sin(c)) -
            std::sin(a) *
                (std::cos(a) * std::pow(std::sin(b), 2) +
                 std::cos(b) * (std::cos(a) * std::cos(b) * std::cos(c) - std::sin(a) * std::sin(c))),
        -(std::sin(a) * (-(std::cos(b) * std::sin(b)) + std::cos(b) * std::cos(c) * std::sin(b))) +
            std::cos(a) * std::sin(b) * std::sin(c),

        std::sin(a) * (std::cos(a) * std::cos(c) - std::cos(b) * std::sin(a) * std::sin(c)) +
            std::cos(a) * (-(std::sin(a) * std::pow(std::sin(b), 2)) +
                           std::cos(b) * (-(std::cos(b) * std::cos(c) * std::sin(a)) -
                                          std::cos(a) * std::sin(c))),
        std::sin(a) * (std::cos(c) * std::sin(a) + std::cos(a) * std::cos(b) * std::sin(c)) +
            std::cos(a) *
                (std::cos(a) * std::pow(std::sin(b), 2) +
                 std::cos(b) * (std::cos(a) * std::cos(b) * std::cos(c) - std::sin(a) * std::sin(c))),
        std::cos(a) * (-(std::cos(b) * std::sin(b)) + std::cos(b) * std::cos(c) * std::sin(b)) +
            std::sin(a) * std::sin(b) * std::sin(c),

        std::cos(b) * std::sin(a) * std::sin(b) +
            std::sin(b) * (-(std::cos(b) * std::cos(c) * std::sin(a)) - std::cos(a) * std::sin(c)),
        -(std::cos(a) * std::cos(b) * std::sin(b)) +
            std::sin(b) * (std::cos(a) * std::cos(b) * std::cos(c) - std::sin(a) * std::sin(c)),
        std::pow(std::cos(b), 2) + std::cos(c) * std::pow(std::sin(b), 2)};
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

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_ROT_A_B_C_HPP
