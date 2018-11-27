// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gablduzz@itp.phys.ethz.ch)
//
// PAPI and time event.

#ifndef DCA_PROFILING_EVENTS_PAPI_AND_TIME_EVENT_HPP
#define DCA_PROFILING_EVENTS_PAPI_AND_TIME_EVENT_HPP

#include <array>
#include <string>
#include <vector>

#include <papi.h>

#include "dca/profiling/events/time_event.hpp"

namespace dca {
namespace profiling {
// dca::profiling::

class PapiAndTimeEvent : public time_event<long long int> {
private:
  constexpr static int max_threads_ = 32;

  constexpr static int nb_time_counter_ = time_event<long long int>::NB_TIME_COUNTERS;
  constexpr static int nb_papi_counter_ = 2;

  // Set of events to record and their name.
  const static std::array<int, nb_papi_counter_> papi_event_types_;
  const static std::array<std::string, nb_papi_counter_> papi_event_names_;

public:
  constexpr static int NB_COUNTERS = nb_time_counter_ + nb_papi_counter_;
  using scalar_type = long long int;

  typedef PapiAndTimeEvent this_type;
  typedef time_event<scalar_type> BaseTimeEvent;

public:
  // The constructor initializes the count.
  PapiAndTimeEvent(std::vector<scalar_type>& counter_ptr_, int id);

  // Read current performance counters and update the total.
  void end();

  // Static initialization.
  static void start();

  // Static cleanup.
  static void stop();

  // Initializes the resources for the current thread.
  static void start_threading(int id);

  // Deallocates the resources for the current thread.
  static void stop_threading(int id);

  // Returns a vector of size NB_COUNTERS with the names of the counters.
  static std::vector<std::string> names();

private:
  static int& papiEventSet(int id);

  std::array<scalar_type, nb_papi_counter_> start_counters_;

  std::vector<scalar_type>* counter_ptr_;

  int thread_id_;
};

inline PapiAndTimeEvent::PapiAndTimeEvent(std::vector<scalar_type>& counter_ref, int id)
    : BaseTimeEvent(counter_ref, id), counter_ptr_(&counter_ref), thread_id_(id) {
  int ret;
  if ((ret = PAPI_read(papiEventSet(thread_id_), start_counters_.data())) != PAPI_OK)
    throw std::logic_error("Error in PAPI_read: " + std::to_string(ret));
}

inline void PapiAndTimeEvent::end() {
  BaseTimeEvent::end();

  std::array<scalar_type, nb_papi_counter_> end_counters;

  int ret;
  if ((ret = PAPI_read(papiEventSet(thread_id_), end_counters.data())) != PAPI_OK) {
      throw std::logic_error("stop: error in PAPI_read: " + std::to_string(ret));
  }

  // Note: the time events stores its counters in positions [0, nb_time_counter_ - 1].
  for (int i = 0; i < nb_papi_counter_; ++i)
    (*counter_ptr_)[nb_time_counter_ + i] += (end_counters[i] - start_counters_[i]);
}

inline int& PapiAndTimeEvent::papiEventSet(int id) {
  static std::vector<int> event_set(max_threads_, PAPI_NULL);
  return event_set[id];
}

}  // profiling
}  // dca

#endif  // DCA_PROFILING_EVENTS_PAPI_AND_TIME_EVENT_HPP
