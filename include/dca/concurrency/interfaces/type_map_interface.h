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

#ifndef DCA_CONCURRENCY_INTERFACES_TYPE_MAP_INTERFACE_H
#define DCA_CONCURRENCY_INTERFACES_TYPE_MAP_INTERFACE_H

#include <complex>
#include "dca/concurrency/concurrency_types.hpp"

namespace dca {
namespace concurrency {
// dca::concurrency::

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY, typename scalar_type>
class type_map_interface {
public:
  static size_t factor() {
    return 1;
  }

  static size_t value() {
    return sizeof(scalar_type);
  }
};

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY, typename scalar_type>
class type_map_interface<LIBRARY, std::complex<scalar_type>> {
public:
  static size_t factor() {
    return 2;
  }

  static size_t value() {
    return sizeof(scalar_type);
  }
};

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_INTERFACES_TYPE_MAP_INTERFACE_H
