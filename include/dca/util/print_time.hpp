// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides functions to print the current time.

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

// Converts the date and time information of 'time_point' to the format dd-mm-yyy hh:mm:ss.
template <typename Clock = std::chrono::system_clock>
std::string print_time(const std::chrono::time_point<Clock>& time_point) {
  std::stringstream s;
  const auto now = Clock::to_time_t(time_point);

  // std::put_time is only available for GCC >= 5.0. We use std:strftime instead.
  // Reference:
  // http://stackoverflow.com/questions/37421747/is-there-a-builtin-alternative-to-stdput-time-for-gcc-5.
  char buffer[24];
  if (0 < std::strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S", std::localtime(&now)))
    s << buffer;

  return s.str();
}

// Returns a string of the current time in the format dd-mm-yyy hh:mm:ss.
template <typename Clock = std::chrono::system_clock>
std::string print_time() {
  return print_time(Clock::now());
}

}  // util
}  // dca

#endif  // DCA_UTIL_PRINT_TIME_HPP
