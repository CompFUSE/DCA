// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gablduzz@itp.phys.ethz.ch)
//
// PAPI and time event.

#include "dca/profiling/events/papi_and_time_event.hpp"

#include <mutex>  // std::call_once.
#include <thread>

namespace dca {
namespace profiling {
// dca::profiling::

const std::array<int, PapiAndTimeEvent::nb_papi_counter_> PapiAndTimeEvent::papi_event_types_{
    PAPI_TOT_INS, PAPI_FP_OPS};
const std::array<std::string, PapiAndTimeEvent::nb_papi_counter_> PapiAndTimeEvent::papi_event_names_{
    "PAPI_TOT_INS", "PAPI_FP_OPS"};

void PapiAndTimeEvent::start() {
  static std::once_flag flag;
  std::call_once(flag, []() {

    int retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT && retval > 0)
      throw std::logic_error("Info :: retval != PAPI_VER_CURRENT && retval>0 .\n");

    if (retval < 0) {
      std::stringstream ss;
      ss << "Info :: retval=" << retval << " .\n";
      throw std::logic_error(ss.str());
    }

    if (PAPI_is_initialized() != PAPI_LOW_LEVEL_INITED)
      throw std::logic_error("Info :: retval != PAPI_LOW_LEVEL_INITED .\n");

    // PAPI_thread_init requires a function that returns an unique id for each thread.
    auto get_id = [] {
      static std::hash<std::thread::id> hash;
      // Map an abstract object of type thread::id to an unique integer value.
      return static_cast<unsigned long>(hash(std::this_thread::get_id()));
    };

    if (PAPI_thread_init(get_id) != PAPI_OK)
      throw std::logic_error("Threaded PAPI initialization failure.");
  });

  start_threading(0);
}

void PapiAndTimeEvent::stop() {
  stop_threading(0);
}

void PapiAndTimeEvent::start_threading(int id) {
  if (id >= max_threads_ || id < 0)
    throw(std::out_of_range("Thread id out of range."));

  PAPI_register_thread();

  if (PAPI_create_eventset(&papiEventSet(id)) != PAPI_OK)
    throw std::logic_error("PAPI_create_eventset");

  for (int i = 0; i < nb_papi_counter_; ++i)
    if (PAPI_add_event(papiEventSet(id), papi_event_types_[i]) != PAPI_OK)
      throw std::logic_error(std::string("PAPI does not find ") + papi_event_names_[i]);

  if (PAPI_start(papiEventSet(id)) != PAPI_OK)
    throw std::logic_error("PAPI does not start");
}

void PapiAndTimeEvent::stop_threading(int id) {
  std::array<scalar_type, nb_papi_counter_> dummy;

  if (PAPI_stop(papiEventSet(id), dummy.data()) != PAPI_OK)
    throw std::logic_error("PAPI does not stop");

  if (PAPI_cleanup_eventset(papiEventSet(id)) != PAPI_OK)
    throw std::logic_error("PAPI_cleanup_eventset");

  if (PAPI_destroy_eventset(&papiEventSet(id)) != PAPI_OK)
    throw std::logic_error("PAPI_destroy_eventset");

  PAPI_unregister_thread();
}

std::vector<std::string> PapiAndTimeEvent::names() {
  std::vector<std::string> names = BaseTimeEvent::names();

  names.insert(names.end(), papi_event_names_.begin(), papi_event_names_.end());
  assert(NB_COUNTERS == names.size());

  return names;
}

}  // profiling
}  // dca
