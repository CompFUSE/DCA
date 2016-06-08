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

#include <cstdlib>
#include "comp_library/profiler_library/profiler.hpp"

class profiling_parameters {
public:
#ifdef DCA_NO_PROFILING
  typedef PROFILER::NullProfiler profiler_type;
#endif  // DCA_NO_PROFILING

#ifdef DCA_COUNTING_PROFILING
  typedef PROFILER::time_event<std::size_t> event_type;
  typedef PROFILER::CountingProfiler<event_type> profiler_type;
#endif  // DCA_COUNTING_PROFILING

#ifdef DCA_PAPI_PROFILING
  typedef PROFILER::papi_and_time_event<long long> event_type;
  typedef PROFILER::CountingProfiler<event_type> profiler_type;
#endif  // DCA_PAPI_PROFILING

#ifdef DCA_CRAYPAT_PROFILING
  typedef std::size_t counter_type;
  typedef dca::CraypatProfiler<counter_type> profiler_type;
#endif  // DCA_CRAYPAT_PROFILING
};

#endif  // PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_PROFILING_PARAMETERS_H
