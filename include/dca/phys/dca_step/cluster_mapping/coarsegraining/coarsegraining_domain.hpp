// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the coarsegraining domain.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_DOMAIN_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_DOMAIN_HPP

#include <stdexcept>
#include <string>
#include <vector>

#include "dca/math/function_transform/domain_specifications.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegrain_domain_names.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
class coarsegraining_domain {
public:
  const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

  typedef double scalar_type;
  typedef std::vector<double> element_type;

  typedef math::transform::domain_specifications<
      scalar_type, element_type, math::transform::DISCRETE, math::transform::KRONECKER_DELTA,
      math::transform::INTERVAL, math::transform::EQUIDISTANT>
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

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_DOMAIN_HPP
