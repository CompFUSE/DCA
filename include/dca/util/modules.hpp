// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
