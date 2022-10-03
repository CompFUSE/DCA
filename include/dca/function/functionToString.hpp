// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file provides a way to do a string dump of a dca::func::function
//

#ifndef DCA_FUNCTION_FUNCTIONTOSTRING_HPP
#define DCA_FUNCTION_FUNCTIONTOSTRING_HPP

#include <string>
#include <sstream>
#include "dca/function/function.hpp"

namespace dca {
namespace func {

template <typename T, class DOMAIN>
void leafTraverse(std::vector<int>& index, const std::vector<std::size_t>& dsizes,
                         const dca::func::function<T, DOMAIN> a_function, const int depth,
                         std::ostringstream& oss) {
  if (depth == 0) {
    oss << '[';
    for (int i = 0; i < dsizes[depth]; ++i) {
      index[0] = i;
      if (i != 0)
        oss << ", ";
      oss << a_function(index);
    }
    oss << "]\n";
  }
  else {
    oss << '[';
    for (int i = 0; i < dsizes[depth]; ++i) {
      index[depth] = i;
      leafTraverse(index, dsizes, a_function, depth - 1, oss);
    }
    oss << "]\n";
  }
}

template <typename T, class DOMAIN>
std::string functionToString(const dca::func::function<T, DOMAIN> a_function) {
  std::ostringstream oss;
  oss << "[";
  auto leaf_sizes = a_function.getDomainSizes();
  std::vector<int> index(leaf_sizes.size());
  leafTraverse(index, leaf_sizes, a_function, leaf_sizes.size() - 1, oss);
  oss << "]\n";
  return oss.str();
}

}  // namespace func
}  // namespace dca
#endif
