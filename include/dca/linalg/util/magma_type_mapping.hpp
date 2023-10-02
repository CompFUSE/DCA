// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W.  Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides better type mapping between host, gpu and magma types
 */

#ifndef DCA_MAGMA_TYPE_MAPPING_HPP
#define DCA_MAGMA_TYPE_MAPPING_HPP

#include <type_traits>
#include <magma_v2.h>
#include "dca/platform/dca_gpu_complex.h"
#include "dca/util/type_mapping.hpp"

namespace dca {
namespace util {

#ifndef DCA_GPU_TYPE_MAPPING_HPP
#define DCA_GPU_TYPE_MAPPING_HPP
template <typename T>
using MAGMAPointerMap = typename std::disjunction<
  OnTypesEqual<T, cuComplex, magmaComplex>,
  OnTypesEqual<T, cuDoubleComplex, magmaDoubleComplex>,
  OnTypesEqual<T, cuComplex*, magmaComplex*>,
  OnTypesEqual<T, cuDoubleComplex*, magmaDoubleComplex*>,
  OnTypesEqual<T, cuComplex**, magmaComplex**>,
  OnTypesEqual<T, cuDoubleComplex**, magmaDoubleComplex**>,
  default_type<void>>::type;

template <typename T>
__device__ __host__ MAGMAPointerMap<T> castMAGMAType(T var) {
  return reinterpret_cast<MAGMAPointerMap<T>>(var);
}

}
}

#endif
