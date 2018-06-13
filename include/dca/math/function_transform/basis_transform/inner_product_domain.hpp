// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the inner product domain used in the basis transformations.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_BASIS_TRANSFORM_INNER_PRODUCT_DOMAIN_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_BASIS_TRANSFORM_INNER_PRODUCT_DOMAIN_HPP

#include <string>
#include <vector>

#include "dca/math/function_transform/boundary_conditions.hpp"
#include "dca/math/function_transform/domain_specifications.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

template <typename dmn_type>
class inner_product_domain {
  typedef typename dmn_type::dmn_specifications_type dmn_specs_type;

public:
  typedef typename dmn_specs_type::scalar_type scalar_type;
  typedef typename dmn_specs_type::element_type element_type;

  const static BOUNDARY_CONDITIONS BOUNDARY_CONDITION = dmn_specs_type::BOUNDARY_CONDITION;

  typedef domain_specifications<scalar_type, element_type, DISCRETE, KRONECKER_DELTA,
                                BOUNDARY_CONDITION, NONEQUIDISTANT>
      dmn_specifications_type;

  typedef inner_product_domain<dmn_type> this_type;

  static bool& is_initialized() {
    static bool initialized = false;
    return initialized;
  }

  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::vector<scalar_type>& get_weights() {
    static std::vector<scalar_type> weights(get_size());
    return weights;
  }

  static std::string& get_name() {
    static std::string name = "inner-product-domain (" + dmn_specifications_type::get_name() + ")";
    return name;
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> elements(get_size());
    return elements;
  }

  static void reset() {
    get_size() = 0;

    get_weights().resize(0);
    get_elements().resize(0);

    is_initialized() = false;
  }

  static void initialize(int level) {
    get_elements().resize(0);

    dmn_type::initialize_integration_domain(level, get_weights(), get_elements());

    get_size() = get_elements().size();
  }
};

}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_BASIS_TRANSFORM_INNER_PRODUCT_DOMAIN_HPP
