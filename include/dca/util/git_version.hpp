// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class stores and prints the most recent git log and current git status.

#ifndef DCA_UTIL_GIT_VERSION_HPP
#define DCA_UTIL_GIT_VERSION_HPP

#include <string>

namespace dca {
namespace util {
// dca::util::

struct GitVersion {
  static const std::string git_log;
  static const std::string git_status;

  static void print();
  static std::string string();
};

}  // util
}  // dca

#endif  // DCA_UTIL_GIT_VERSION_HPP
