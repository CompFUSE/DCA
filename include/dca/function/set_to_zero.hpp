// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Bart Ydens
//
// This class assures that all functions with a scalartype of float, double, std::complex<float>, or
// std::complex<double> are initialized to zero.
// Functions with a different 'scalartype' are not initialized.
//
// TODO: Rewrite/remove with rewriting of function class.

#ifndef DCA_FUNCTION_SET_TO_ZERO_HPP
#define DCA_FUNCTION_SET_TO_ZERO_HPP

#include <complex>

namespace dca {
namespace func {
// dca::func::

struct set_to_zero {
  template <class whatever_t>
  static void execute(whatever_t& whatever);
};

template <class whatever_t>
inline void set_to_zero::execute(whatever_t& /*whatever*/) {}

template <>
inline void set_to_zero::execute(short& whatever) {
  whatever = 0;
}

template <>
inline void set_to_zero::execute(int& whatever) {
  whatever = 0;
}

template <>
inline void set_to_zero::execute(long& whatever) {
  whatever = 0;
}

template <>
inline void set_to_zero::execute(size_t& whatever) {
  whatever = 0;
}

template <>
inline void set_to_zero::execute(float& whatever) {
  whatever = 0.;
}

template <>
inline void set_to_zero::execute(double& whatever) {
  whatever = 0.;
}

template <>
inline void set_to_zero::execute(std::complex<float>& whatever) {
  whatever = 0.;
}

template <>
inline void set_to_zero::execute(std::complex<double>& whatever) {
  whatever = 0.;
}

}  // func
}  // dca

#endif  // DCA_FUNCTION_SET_TO_ZERO_HPP
