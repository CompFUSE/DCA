// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class stores and prints a list of the currently loaded modules.

#ifndef DCA_UTIL_GIT_MODULES_HPP
#define DCA_UTIL_GIT_MODULES_HPP

#include <string>

namespace dca {
namespace util {
// dca::util::

struct Modules {
  static const std::string module_list;

  static void print();
  static std::string string();
};

}  // util
}  // dca

#endif  // DCA_UTIL_GIT_MODULES_HPP
