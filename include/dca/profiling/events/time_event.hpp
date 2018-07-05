// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Time event.

#ifndef DCA_PROFILING_EVENTS_TIME_EVENT_HPP
#define DCA_PROFILING_EVENTS_TIME_EVENT_HPP

#include <cassert>
#include <string>
#include <vector>

#include "dca/profiling/events/time.hpp"

namespace dca {
namespace profiling {
// dca::profiling::

template <typename scalartype>
class time_event {
public:
  typedef scalartype scalar_type;

  const static int NB_TIME_COUNTERS = 8;

  const static int NB_COUNTERS = NB_TIME_COUNTERS;

  time_event(std::vector<scalartype>& my_counter, int id);

  void end();

  static void start() {}
  static void stop() {}

  static void start_threading(int /*id*/) {}
  static void stop_threading(int /*id*/) {}

  static std::vector<std::string> names();

  static void normalize(std::vector<scalartype>& counters);

  void update_counter(std::vector<scalartype>& counters, Duration wallDuration,
                      Duration userDuration, Duration systemDuration);

private:
  WallTime wallTime;
  UserTime userTime;
  SystemTime systemTime;

  std::vector<scalartype>& ourCounter;

  int thread_id;
};

template <typename scalartype>
time_event<scalartype>::time_event(std::vector<scalartype>& ourCountr, int id)
    : wallTime(), userTime(), systemTime(userTime.usage), ourCounter(ourCountr), thread_id(id) {}

template <typename scalartype>
void time_event<scalartype>::end() {
  WallTime endWallTime;
  UserTime endUserTime;
  SystemTime endSystemTime(endUserTime.usage);

  Duration wallDuration(endWallTime, wallTime);
  Duration userDuration(endUserTime, userTime);
  Duration systemDuration(endSystemTime, systemTime);

  update_counter(ourCounter, wallDuration, userDuration, systemDuration);
}

template <typename scalartype>
void time_event<scalartype>::update_counter(std::vector<scalartype>& counters, Duration wallDuration,
                                            Duration userDuration, Duration systemDuration) {
  Duration previousWallDuration(counters[0], counters[1]);
  Duration previousUserDuration(counters[2], counters[3]);
  Duration previousSystemDuration(counters[4], counters[5]);

  Duration newWallDuration(previousWallDuration, wallDuration);
  Duration newUserDuration(previousUserDuration, userDuration);
  Duration newSystemDuration(previousSystemDuration, systemDuration);

  counters[0] = newWallDuration.sec;
  counters[1] = newWallDuration.usec;

  counters[2] = newUserDuration.sec;
  counters[3] = newUserDuration.usec;

  counters[4] = newSystemDuration.sec;
  counters[5] = newSystemDuration.usec;

  counters[6] += 1;  // number of calls
  counters[7] = 1;   // number of threads
}

template <typename scalartype>
void time_event<scalartype>::normalize(std::vector<scalartype>& counters) {
  Duration wallDuration(counters[0], counters[1]);
  Duration userDuration(counters[2], counters[3]);
  Duration systemDuration(counters[4], counters[5]);

  counters[0] = wallDuration.sec;
  counters[1] = wallDuration.usec;

  counters[2] = userDuration.sec;
  counters[3] = userDuration.usec;

  counters[4] = systemDuration.sec;
  counters[5] = systemDuration.usec;
}

template <typename scalartype>
std::vector<std::string> time_event<scalartype>::names() {
  std::vector<std::string> names(0);

  names.push_back("wall        seconds");  // 0
  names.push_back("wall   microseconds");  // 1
  names.push_back("user        seconds");  // 2
  names.push_back("user   microseconds");  // 3
  names.push_back("system      seconds");  // 4
  names.push_back("system microseconds");  // 5
  names.push_back("calls              ");  // 6
  names.push_back("threads            ");  // 7

  assert(NB_COUNTERS == names.size());

  return names;
}

}  // profiling
}  // dca

#endif  // DCA_PROFILING_EVENTS_TIME_EVENT_HPP
