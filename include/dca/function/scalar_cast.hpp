// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides utility functions that cast real/complex numbers to real/complex numbers of a
// potentially different floating point type.
//
// TODO: These functions are only used in function::slice and function::distribute. One should check
//       whether the cast complex to real is actually carried out. If not, the remaining functions
//       could be replaced with inline static_cast's.

#ifndef DCA_FUNCTION_SCALAR_CAST_HPP
#define DCA_FUNCTION_SCALAR_CAST_HPP

#include <complex>

namespace dca {
namespace func {
// dca::func::

// Real/complex --> real
template <typename ScalarType>
struct ScalarCast {
  template <typename OtherScalarType>
  static ScalarType execute(OtherScalarType value) {
    return static_cast<ScalarType>(value);
  }

  template <typename OtherScalarType>
  static ScalarType execute(std::complex<OtherScalarType> value) {
    return static_cast<ScalarType>(value.real());
  }
};

// Real/complex --> complex
template <typename ScalarType>
struct ScalarCast<std::complex<ScalarType>> {
  template <typename OtherScalarType>
  static std::complex<ScalarType> execute(OtherScalarType value) {
    return static_cast<std::complex<ScalarType>>(value);
  }

  template <typename OtherScalarType>
  static std::complex<ScalarType> execute(std::complex<OtherScalarType> value) {
    return static_cast<std::complex<ScalarType>>(value);
  }
};

}  // func
}  // dca

#endif  // DCA_FUNCTION_SCALAR_CAST_HPP
