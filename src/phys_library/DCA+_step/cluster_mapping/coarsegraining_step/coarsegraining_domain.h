// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_DOMAIN_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_DOMAIN_H

#include <stdexcept>
#include <string>
#include <vector>

#include "math_library/functional_transforms/domain_specifications/domain_specifications.hpp"
#include "phys_library/DCA+_step/cluster_mapping/coarsegraining_step/coarsegraining_names.hpp"

namespace DCA {
template <typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
class coarsegraining_domain {
public:
  const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

  typedef double scalar_type;
  typedef std::vector<double> element_type;

  typedef math_algorithms::domain_specifications<
      scalar_type, element_type, math_algorithms::DISCRETE, math_algorithms::KRONECKER_DELTA,
      math_algorithms::INTERVAL, math_algorithms::EQUIDISTANT>
      dmn_specifications_type;

public:
  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::string get_name() {
    static std::string name = "coarsegrain_domain (" + to_str(NAME) + ")";
    return name;
  }

  static std::vector<scalar_type>& get_weights() {
    static std::vector<scalar_type> weights(0);
    return weights;
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> elements(0);
    return elements;
  }

  static void set_elements(int K_ind) {
    switch (NAME) {
      case K: {
        std::vector<element_type> elements = coarsegraining_domain<K_dmn, ORIGIN>::get_elements();

        for (int q_ind = 0; q_ind < elements.size(); q_ind++)
          for (int d_ind = 0; d_ind < DIMENSION; d_ind++)
            elements[q_ind][d_ind] += K_dmn::get_elements()[K_ind][d_ind];

        get_elements() = elements;
      } break;

      case TETRAHEDRON_K: {
        std::vector<element_type> elements =
            coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_elements();

        for (int q_ind = 0; q_ind < elements.size(); q_ind++)
          for (int d_ind = 0; d_ind < DIMENSION; d_ind++)
            elements[q_ind][d_ind] += K_dmn::get_elements()[K_ind][d_ind];

        get_elements() = elements;
      } break;

      default:
        throw std::logic_error(__FUNCTION__);
    }
  }

  static void set_elements(int K_ind, int Q_ind) {
    std::vector<element_type> elements = coarsegraining_domain<K_dmn, ORIGIN>::get_elements();

    switch (NAME) {
      case K_PLUS_Q: {
        for (int q_ind = 0; q_ind < elements.size(); q_ind++) {
          for (int d_ind = 0; d_ind < DIMENSION; d_ind++) {
            elements[q_ind][d_ind] += K_dmn::get_elements()[K_ind][d_ind];
            elements[q_ind][d_ind] += K_dmn::get_elements()[Q_ind][d_ind];
          }
        }
      } break;

      case Q_MINUS_K: {
        for (int q_ind = 0; q_ind < elements.size(); q_ind++) {
          for (int d_ind = 0; d_ind < DIMENSION; d_ind++) {
            elements[q_ind][d_ind] *= -1;
            elements[q_ind][d_ind] -= K_dmn::get_elements()[K_ind][d_ind];
            elements[q_ind][d_ind] += K_dmn::get_elements()[Q_ind][d_ind];
          }
        }
      } break;

      default:
        throw std::logic_error(__FUNCTION__);
    }

    get_elements() = elements;
  }
};
}

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_STEP_COARSEGRAINING_DOMAIN_H
