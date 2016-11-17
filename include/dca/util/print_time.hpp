// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides a function that prints the current time.

#ifndef DCA_UTIL_PRINT_TIME_HPP
#define DCA_UTIL_PRINT_TIME_HPP

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

namespace dca {
namespace util {
// dca::util::

template <typename clock = std::chrono::system_clock>
std::string print_time() {
  std::stringstream s;
  auto now = clock::to_time_t(clock::now());

  // std::put_time is only available for GCC >= 5.0. We use std:strftime instead.
  // Reference:
  // http://stackoverflow.com/questions/37421747/is-there-a-builtin-alternative-to-stdput-time-for-gcc-5.
  char buffer[24];
  if (0 < std::strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S", std::localtime(&now)))
    s << buffer;

  return s.str();
}

}  // util
}  // dca

#endif  // DCA_UTIL_PRINT_TIME_HPP
