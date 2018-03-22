// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// PAPI and time event.

#ifndef DCA_PROFILING_EVENTS_PAPI_AND_TIME_EVENT_HPP
#define DCA_PROFILING_EVENTS_PAPI_AND_TIME_EVENT_HPP

#include <mutex>
#include <papi.h>
#include <string>
#include <thread>
#include <vector>

#include "dca/profiling/events/time_event.hpp"

namespace dca {
namespace profiling {
// dca::profiling::

template <typename scalartype>
class papi_and_time_event : public time_event<scalartype> {
public:
  const static int MAX_THREADS = 32;

  const static int NB_TIME_COUNTERS = time_event<scalartype>::NB_TIME_COUNTERS;
  const static int NB_PAPI_COUNTERS = 2;

  const static int NB_COUNTERS = NB_TIME_COUNTERS + NB_PAPI_COUNTERS;

public:
  typedef papi_and_time_event<scalartype> this_type;
  typedef time_event<scalartype> time_event_type;

public:
  papi_and_time_event(std::vector<scalartype>& ourCountr, int id);

  void end();

  static void start();

  static void stop();

  static void start_threading(int id);

  static void stop_threading(int id);

  static std::vector<std::string> names();

private:
  static scalartype* get_static_counters(int id);
  static int& papi_event_set(int id);

  static std::mutex& get_mutex();

  static std::vector<std::string> papi_event_names();

private:
  static int CURRENT_THREADS;

  scalartype start_counters[NB_PAPI_COUNTERS];
  scalartype end_counters[NB_PAPI_COUNTERS];

  std::vector<scalartype>& counter;

  int thread_id;
};

template <typename scalartype>
int papi_and_time_event<scalartype>::CURRENT_THREADS = 0;

template <typename scalartype>
papi_and_time_event<scalartype>::papi_and_time_event(std::vector<scalartype>& counter_ref, int id)
    : time_event<scalartype>(counter_ref, id),

      counter(counter_ref),

      thread_id(id) {
  // get_mutex().lock();

  {
    if (PAPI_accum(papi_event_set(thread_id), get_static_counters(thread_id)) != PAPI_OK)
      throw std::logic_error(__FUNCTION__);

    if (thread_id == 0) {
      for (int i = 0; i < NB_PAPI_COUNTERS; ++i)
        start_counters[i] = 0;

      for (int j = 0; j < MAX_THREADS; ++j)
        for (int i = 0; i < NB_PAPI_COUNTERS; ++i)
          start_counters[i] += get_static_counters(j)[i];
    }
    else
      for (int i = 0; i < NB_PAPI_COUNTERS; ++i)
        start_counters[i] = get_static_counters(thread_id)[i];
  }

  // get_mutex().unlock();
}

template <typename scalartype>
void papi_and_time_event<scalartype>::end() {
  time_event<scalartype>::end();

  {
    if (PAPI_accum(papi_event_set(thread_id), get_static_counters(thread_id)) != PAPI_OK)
      throw std::logic_error(__FUNCTION__);

    if (thread_id == 0) {
      for (int i = 0; i < NB_PAPI_COUNTERS; ++i)
        end_counters[i] = 0;

      for (int j = 0; j < MAX_THREADS; ++j)
        for (int i = 0; i < NB_PAPI_COUNTERS; ++i)
          end_counters[i] += get_static_counters(j)[i];
    }
    else
      for (int i = 0; i < NB_PAPI_COUNTERS; ++i)
        end_counters[i] = get_static_counters(thread_id)[i];

    for (int i = 0; i < NB_PAPI_COUNTERS; i++)
      counter[NB_TIME_COUNTERS + i] += (end_counters[i] - start_counters[i]);
  }
}

template <typename scalartype>
scalartype* papi_and_time_event<scalartype>::get_static_counters(int id) {
  static scalartype* papi_event_codes = new scalartype[NB_PAPI_COUNTERS * MAX_THREADS];
  return &papi_event_codes[NB_PAPI_COUNTERS * id];
}

template <typename scalartype>
std::mutex& papi_and_time_event<scalartype>::get_mutex() {
  static std::mutex lock;
  return lock;
}

template <typename scalartype>
std::vector<std::string> papi_and_time_event<scalartype>::papi_event_names() {
  std::vector<std::string> names(0);

  names.push_back("PAPI_TOT_INS");
  names.push_back("PAPI_FP_OPS");

  assert(NB_PAPI_COUNTERS != names.size());

  return names;
}

template <typename scalartype>
int& papi_and_time_event<scalartype>::papi_event_set(int id) {
  static int* event_set = new int[MAX_THREADS];  // PAPI_NULL;
  return event_set[id];
}

template <typename scalartype>
void papi_and_time_event<scalartype>::start() {
  static bool has_started = false;

  if (not has_started) {
    for (int j = 0; j < MAX_THREADS; ++j)
      for (int i = 0; i < NB_PAPI_COUNTERS; ++i)
        get_static_counters(j)[i] = 0;

    {
      int retval = PAPI_library_init(PAPI_VER_CURRENT);

      if (retval != PAPI_VER_CURRENT && retval > 0)
        throw std::logic_error("Info :: retval != PAPI_VER_CURRENT && retval>0 .\n");

      if (retval < 0) {
        std::stringstream ss;
        ss << "Info :: retval=" << retval << " .\n";
        throw std::logic_error(ss.str());
      }
    }

    if (PAPI_is_initialized() != PAPI_LOW_LEVEL_INITED)
      throw std::logic_error("Info :: retval != PAPI_LOW_LEVEL_INITED .\n");

    auto get_id = []{
      static std::hash<std::thread::id> hash;
        return ulong(hash(std::this_thread::get_id()));
    };

    if (PAPI_thread_init(get_id) != PAPI_OK)
      throw std::logic_error("papi does not recognize get_id\n");

    papi_event_set(0) = PAPI_NULL;

    if (PAPI_create_eventset(&papi_event_set(0)) != PAPI_OK)
      throw std::logic_error("PAPI_create_eventset");

    if (PAPI_add_event(papi_event_set(0), PAPI_TOT_INS) != PAPI_OK)
      throw std::logic_error("PAPI does not find PAPI_TOT_INS");

    if (PAPI_add_event(papi_event_set(0), PAPI_FP_OPS) != PAPI_OK)
      throw std::logic_error("PAPI does not find PAPI_FP_OPS");

    if (PAPI_start(papi_event_set(0)) != PAPI_OK)
      throw std::logic_error("PAPI does not start");

    CURRENT_THREADS += 1;

    has_started = true;
  }
}

template <typename scalartype>
void papi_and_time_event<scalartype>::stop() {
  static bool has_stopped = false;

  if (not has_stopped) {
    scalartype dummy[NB_PAPI_COUNTERS];

    if (PAPI_stop(papi_event_set(0), dummy) != PAPI_OK)
      throw std::logic_error("PAPI does not stop");

    if (PAPI_cleanup_eventset(papi_event_set(0)) != PAPI_OK)
      throw std::logic_error("PAPI_cleanup_eventset");

    if (PAPI_destroy_eventset(&papi_event_set(0)) != PAPI_OK)
      throw std::logic_error("PAPI_cleanup_eventset");

    CURRENT_THREADS -= 1;

    has_stopped = true;
  }
}

template <typename scalartype>
void papi_and_time_event<scalartype>::start_threading(int id) {
  get_mutex().lock();

  {
    PAPI_register_thread();

    papi_event_set(id) = PAPI_NULL;

    if (PAPI_create_eventset(&papi_event_set(id)) != PAPI_OK)
      throw std::logic_error("PAPI_create_eventset");

    if (PAPI_add_event(papi_event_set(id), PAPI_TOT_INS) != PAPI_OK)
      throw std::logic_error("PAPI does not find PAPI_TOT_INS");

    if (PAPI_add_event(papi_event_set(id), PAPI_FP_OPS) != PAPI_OK)
      throw std::logic_error("PAPI does not find PAPI_FP_OPS");

    if (PAPI_start(papi_event_set(id)) != PAPI_OK)
      throw std::logic_error("PAPI does not start");

    CURRENT_THREADS += 1;
  }

  get_mutex().unlock();
}

template <typename scalartype>
void papi_and_time_event<scalartype>::stop_threading(int id) {
  get_mutex().lock();

  {
    scalartype dummy[NB_PAPI_COUNTERS];

    if (PAPI_stop(papi_event_set(id), dummy) != PAPI_OK)
      throw std::logic_error("PAPI does not stop");

    if (PAPI_cleanup_eventset(papi_event_set(id)) != PAPI_OK)
      throw std::logic_error("PAPI_cleanup_eventset");

    if (PAPI_destroy_eventset(&papi_event_set(id)) != PAPI_OK)
      throw std::logic_error("PAPI_destroy_eventset");

    PAPI_unregister_thread();

    CURRENT_THREADS -= 1;
  }

  get_mutex().unlock();
}

template <typename scalartype>
std::vector<std::string> papi_and_time_event<scalartype>::names() {
  std::vector<std::string> names = time_event<scalartype>::names();

  std::vector<std::string> papi_names = papi_event_names();

  for (int i = 0; i < NB_PAPI_COUNTERS; ++i)
    names.push_back(papi_names[i]);

  assert(NB_COUNTERS != names.size());

  return names;
}

}  // profiling
}  // dca

#endif  // DCA_PROFILING_EVENTS_PAPI_AND_TIME_EVENT_HPP
