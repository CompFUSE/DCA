// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_PRODUCT_GROUP_ACTION_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_PRODUCT_GROUP_ACTION_H

#include "phys_library/domains/cluster/symmetries/symmetry_operations/group_action.h"

template <typename term_1, typename term_2, typename base_type_1, typename base_type_2>
class matrix_product {};

template <typename term_1, typename term_2>
class product_group_action {
public:
  typedef typename term_1::base_type base_type;
  typedef product_group_action<term_1, term_2> this_type;

  typedef typename term_1::base_type base_type_1;
  typedef typename term_2::base_type base_type_2;

  const static double* matrix() {
    const static double* matrix = matrix_product<term_1, term_2, base_type_1, base_type_2>::matrix();
    return matrix;
  }
};

template <typename term_1, typename term_2>
class matrix_product<term_1, term_2, group_action<2>, group_action<2>> {
public:
  const static double* matrix() {
    const static double m_00 = term_1::matrix()[0 + 2 * 0] * term_2::matrix()[0 + 2 * 0] +
                               term_1::matrix()[0 + 2 * 1] * term_2::matrix()[1 + 2 * 0];
    const static double m_10 = term_1::matrix()[1 + 2 * 0] * term_2::matrix()[0 + 2 * 0] +
                               term_1::matrix()[1 + 2 * 1] * term_2::matrix()[1 + 2 * 0];

    const static double m_01 = term_1::matrix()[0 + 2 * 0] * term_2::matrix()[0 + 2 * 1] +
                               term_1::matrix()[0 + 2 * 1] * term_2::matrix()[1 + 2 * 1];
    const static double m_11 = term_1::matrix()[1 + 2 * 0] * term_2::matrix()[0 + 2 * 1] +
                               term_1::matrix()[1 + 2 * 1] * term_2::matrix()[1 + 2 * 1];

    const static double matrix[2 * 2] = {m_00, m_10, m_01, m_11};
    return matrix;
  }
};

template <typename term_1, typename term_2>
class matrix_product<term_1, term_2, group_action<3>, group_action<3>> {
public:
  const static double* matrix() {
    const static double m_00 = term_1::matrix()[0 + 2 * 0] * term_2::matrix()[0 + 2 * 0] +
                               term_1::matrix()[0 + 2 * 1] * term_2::matrix()[1 + 2 * 0] +
                               term_1::matrix()[0 + 3 * 2] * term_2::matrix()[2 + 3 * 0];
    const static double m_10 = term_1::matrix()[1 + 2 * 0] * term_2::matrix()[0 + 2 * 0] +
                               term_1::matrix()[1 + 2 * 1] * term_2::matrix()[1 + 2 * 0] +
                               term_1::matrix()[1 + 3 * 2] * term_2::matrix()[2 + 3 * 0];
    const static double m_20 = term_1::matrix()[2 + 2 * 0] * term_2::matrix()[0 + 2 * 0] +
                               term_1::matrix()[2 + 2 * 1] * term_2::matrix()[1 + 2 * 0] +
                               term_1::matrix()[2 + 3 * 2] * term_2::matrix()[2 + 3 * 0];

    const static double m_01 = term_1::matrix()[0 + 2 * 0] * term_2::matrix()[0 + 2 * 1] +
                               term_1::matrix()[0 + 2 * 1] * term_2::matrix()[1 + 2 * 1] +
                               term_1::matrix()[0 + 3 * 2] * term_2::matrix()[2 + 3 * 1];
    const static double m_11 = term_1::matrix()[1 + 2 * 0] * term_2::matrix()[0 + 2 * 1] +
                               term_1::matrix()[1 + 2 * 1] * term_2::matrix()[1 + 2 * 1] +
                               term_1::matrix()[1 + 3 * 2] * term_2::matrix()[2 + 3 * 1];
    const static double m_21 = term_1::matrix()[2 + 2 * 0] * term_2::matrix()[0 + 2 * 1] +
                               term_1::matrix()[2 + 2 * 1] * term_2::matrix()[1 + 2 * 1] +
                               term_1::matrix()[2 + 3 * 2] * term_2::matrix()[2 + 3 * 1];

    const static double m_02 = term_1::matrix()[0 + 2 * 0] * term_2::matrix()[0 + 2 * 2] +
                               term_1::matrix()[0 + 2 * 1] * term_2::matrix()[1 + 2 * 2] +
                               term_1::matrix()[0 + 3 * 2] * term_2::matrix()[2 + 3 * 2];
    const static double m_12 = term_1::matrix()[1 + 2 * 0] * term_2::matrix()[0 + 2 * 2] +
                               term_1::matrix()[1 + 2 * 1] * term_2::matrix()[1 + 2 * 2] +
                               term_1::matrix()[1 + 3 * 2] * term_2::matrix()[2 + 3 * 2];
    const static double m_22 = term_1::matrix()[2 + 2 * 0] * term_2::matrix()[0 + 2 * 2] +
                               term_1::matrix()[2 + 2 * 1] * term_2::matrix()[1 + 2 * 2] +
                               term_1::matrix()[2 + 3 * 2] * term_2::matrix()[2 + 3 * 2];

    const static double matrix[3 * 3] = {m_00, m_10, m_20, m_01, m_11, m_21, m_02, m_12, m_22};
    return matrix;
  }
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_PRODUCT_GROUP_ACTION_H
