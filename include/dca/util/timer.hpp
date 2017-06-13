// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// The class allows to easily time methods or blocks.

#ifndef DCA_UTIL_TIMER_HPP
#define DCA_UTIL_TIMER_HPP

#include <chrono>
#include <string>

namespace dca {
namespace util {
// dca::util::

class Timer {
public:
  // Prints the current time.
  Timer(const std::string& name);

  // Prints the current time and the time passed since construction.
  ~Timer();

private:
  const std::string name_;
  const std::chrono::time_point<std::chrono::system_clock> start_;
};

}  // util
}  // dca

#endif  // DCA_UTIL_TIMER_HPP
