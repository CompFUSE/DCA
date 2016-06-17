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
#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_INNER_PRODUCT_DOMAIN_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_INNER_PRODUCT_DOMAIN_H

#include <string>
#include <vector>

#include "math_library/functional_transforms/domain_specifications/domain_specifications.hpp"
#include "math_library/functional_transforms/typedefs.hpp"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::

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

public:
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

}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_TRANSFORMS_INNER_PRODUCT_DOMAIN_H
