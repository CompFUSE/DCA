// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements a simple, but generic domain that can be used with the function class.
// It's templated on its size and the element type. If needed, elements can be assigned and stored.
//
// TODO: The size type should be std::size_t. However, since the underlying type of std::size_t is
//       compiler dependent, the name of the domain would be compiler dependent, too. This would
//       make writing tests a little harder.

#ifndef DCA_FUNCTION_DOMAINS_DMN_HPP
#define DCA_FUNCTION_DOMAINS_DMN_HPP

#include <stdexcept>
#include <string>
#include <vector>

#include "dca/util/type_utils.hpp"  // for dca::util::type_name

namespace dca {
namespace func {
// dca::func::

template <int size, class element_t = int>
class dmn {
public:
  static_assert(size > 0, "Size must pe positive.");

  using element_type = element_t;
  using this_type = dmn<size, element_t>;

  static int get_size() {
    return size;
  }
  static int dmn_size() {
    return get_size();
  }

  static std::string get_name() {
    return dca::util::type_name<this_type>();
  }

  static void set_elements(const std::vector<element_t>& elements) {
    if (elements.size() != size)
      throw std::logic_error(
          "The size of the passed elements doesn't match the size of the domain.");
    elements_ = elements;
  }

  // For the moment cannot return a const reference since dmn_0's get_elements only returns a
  // non-const reference.
  static std::vector<element_t>& get_elements() {
    if (elements_.size() != size)
      throw std::logic_error("Elements have not been set.");
    return elements_;
  }

private:
  static std::vector<element_t> elements_;
};

template <int size, class element_t>
std::vector<element_t> dmn<size, element_t>::elements_;

}  // func
}  // dca

#endif  // DCA_FUNCTION_DOMAINS_DMN_HPP
