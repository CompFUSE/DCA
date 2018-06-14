// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Mirror alpha beta.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_MIRR_A_B_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_MIRR_A_B_HPP

#include <cmath>
#include "dca/phys/domains/cluster/symmetries/symmetry_operations/group_action.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <int alpha_num, int alpha_den, int beta_num, int beta_den>
class mirr_a_b : public group_action<3> {
public:
  typedef group_action<3> base_type;
  typedef mirr_a_b<alpha_num, alpha_den, beta_num, beta_den> this_type;

  mirr_a_b(){};

  const static double* matrix() {
    double a = 2 * M_PI * double(alpha_num) / double(alpha_den);  // polar coordinate \theta
    double b = 2 * M_PI * double(beta_num) / double(beta_den);    // polar coordinate \phi

    // matrix in column-major format...
    static double matrix[3 * 3] = {
        std::pow(std::cos(a), 2) -
            std::sin(a) *
                (-(std::pow(std::cos(b), 2) * std::sin(a)) + std::sin(a) * std::pow(std::sin(b), 2)),
        std::cos(a) * std::sin(a) -
            std::sin(a) *
                (std::cos(a) * std::pow(std::cos(b), 2) - std::cos(a) * std::pow(std::sin(b), 2)),
        -2 * std::cos(b) * std::sin(a) * std::sin(b),

        std::cos(a) * std::sin(a) +
            std::cos(a) *
                (-(std::pow(std::cos(b), 2) * std::sin(a)) + std::sin(a) * std::pow(std::sin(b), 2)),
        std::pow(std::sin(a), 2) +
            std::cos(a) *
                (std::cos(a) * std::pow(std::cos(b), 2) - std::cos(a) * std::pow(std::sin(b), 2)),
        2 * std::cos(a) * std::cos(b) * std::sin(b),

        -2 * std::cos(b) * std::sin(a) * std::sin(b), 2 * std::cos(a) * std::cos(b) * std::sin(b),
        -std::pow(std::cos(b), 2) + std::pow(std::sin(b), 2)};
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
  // P := ({
  //    {1, 0, 0},
  //    {0, 1, 0},
  //    {0, 0, -1}
  //   })
  //
  // Transpose[
  //   Z[\[Alpha]].X[\[Beta]].P.X[-\[Beta]].Z[-\[Alpha]]]
  // matrix = CForm[Transpose[
  //    X[\[Alpha]].Z[\[Beta]].P.Z[-\[Beta]].X[-\[Alpha]]] /. {\
  // \[Alpha] -> a, \[Beta] -> b, }]
  */
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_MIRR_A_B_HPP
