// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// This class provides type defintions for profiling.

#ifndef PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_PROFILING_PARAMETERS_H
#define PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_PROFILING_PARAMETERS_H

#ifdef NO_PROFILING
#include "comp_library/profiler_library/profilers/null_profiler.h"
#endif
#if defined(COUNTING_PROFILING) || defined(DCA_PAPI_PROFILING)
#include "comp_library/profiler_library/profilers/counting_profiler.h"
#endif

template <class concurrency_type>
class profiling_parameters {
public:
#ifdef NO_PROFILING
  typedef PROFILER::no_profiling_mode profiler_type;
#endif

#ifdef COUNTING_PROFILING
  typedef PROFILER::time_event<size_t> event_type;
  typedef PROFILER::CountingProfiler<event_type> profiler_type;
#endif

#ifdef DCA_CRAYPAT_PROFILING
  typedef size_t counter_type;
  typedef dca::CraypatProfiler<counter_type> profiler_type;
#endif

#ifdef DCA_PAPI_PROFILING
  typedef PROFILER::papi_and_time_event<long long> event_type;
  typedef PROFILER::CountingProfiler<event_type> profiler_type;
#endif
};

#endif  // PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_PROFILING_PARAMETERS_H
