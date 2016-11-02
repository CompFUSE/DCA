// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Elements of the dihedral groups D_n(m).

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_SN_3D_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_SN_3D_HPP

#include <cmath>
#include "dca/phys/domains/cluster/symmetries/symmetry_operations/group_action.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <int axis, int n, int m>
class Sn_3D : public group_action<3> {};

template <int n, int m>
class Sn_3D<0, n, m> : public group_action<3> {
public:
  typedef group_action<3> base_type;
  typedef Sn_3D<0, n, m> this_type;

  Sn_3D(){};

  const static double* matrix() {
    double c(std::cos(2 * M_PI * double(n) / double(m)));
    double s(std::sin(2 * M_PI * double(n) / double(m)));

    static double matrix[3 * 3] = {
        1., 0., 0., 0., c * c - s * s, 2 * s * c, 0., 2 * s * c, -(c * c - s * s)};

    return matrix;
  }
};

template <int n, int m>
class Sn_3D<1, n, m> : public group_action<3> {
public:
  typedef group_action<3> base_type;
  typedef Sn_3D<1, n, m> this_type;

  Sn_3D(){};

  const static double* matrix() {
    double c(std::cos(2 * M_PI * double(n) / double(m)));
    double s(std::sin(2 * M_PI * double(n) / double(m)));

    static double matrix[3 * 3] = {c * c - s * s, 0., 2 * s * c,       0., 1., 0.,
                                   2 * s * c,     0., -(c * c - s * s)};

    return matrix;
  }
};

template <int n, int m>
class Sn_3D<2, n, m> : public group_action<3> {
public:
  typedef group_action<3> base_type;
  typedef Sn_3D<2, n, m> this_type;

  Sn_3D(){};

  const static double* matrix() {
    double c(std::cos(2 * M_PI * double(n) / double(m)));
    double s(std::sin(2 * M_PI * double(n) / double(m)));

    static double matrix[3 * 3] = {
        c * c - s * s, 2 * s * c, 0., 2 * s * c, -(c * c - s * s), 0., 0., 0., 1.};

    return matrix;
  }
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_SN_3D_HPP
