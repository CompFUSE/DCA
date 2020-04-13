// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file provides a utility method to compute the difference of two functions.

#ifndef DCA_FUNCTION_UTIL_DIFFERENCE_HPP
#define DCA_FUNCTION_UTIL_DIFFERENCE_HPP

#include <cassert>
#include <cmath>
#include <type_traits>

#include "dca/function/function.hpp"

namespace dca {
namespace func {
namespace util {
// dca::func::util::

struct Difference {
  double l1, l2, l_inf;
};

// Returns the l1, l2, and l_inf relative (w.r.t. f1) differences of f1 and f2.
template <typename Scalartype1, typename Scalartype2, class Dmn1, class Dmn2>
Difference difference(const function<Scalartype1, Dmn1>& f1, const function<Scalartype2, Dmn2>& f2) {
  if (f1.size() != f2.size())
    throw(std::logic_error("The two functions have different size."));

  double l1 = 0.;
  double l2 = 0.;
  double linf = 0.;

  double l1_error = 0.;
  double l2_error = 0.;
  double linf_error = 0.;

  using DiffType = std::common_type_t<Scalartype1, Scalartype2>;

  for (int i = 0; i < f2.size(); ++i) {
    const double f_abs = std::abs(f1(i));
    l1 += f_abs;
    l2 += f_abs * f_abs;
    linf = std::max(linf, f_abs);

    const auto diff = static_cast<DiffType>(f1(i)) - static_cast<DiffType>(f2(i));
    const double err = std::abs(diff);
    l1_error += err;
    l2_error += err * err;
    linf_error = std::max(linf_error, err);
  }

  l1_error /= l1;
  l2_error = std::sqrt(l2_error / l2);
  linf_error /= linf;

  return Difference{l1_error, l2_error, linf_error};
}

}  // util
}  // func
}  // dca

#endif  // DCA_FUNCTION_UTIL_DIFFERENCE_HPP
