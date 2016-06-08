// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// This file provides a function that prints the current time.
//
// TODO: Make this file a source file?

#ifndef DCA_UTIL_PRINT_TIME_HPP
#define DCA_UTIL_PRINT_TIME_HPP

#include <string>

namespace dca {
namespace util {
// dca::util::

std::string print_time() {
  time_t rawtime;
  struct tm* timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);

  char buffer[80];
  strftime(buffer, 80, "%d-%m-%Y %I:%M:%S", timeinfo);
  std::string str(buffer);

  return str;
}

}  // util
}  // dca

#endif  // DCA_UTIL_PRINT_TIME_HPP
