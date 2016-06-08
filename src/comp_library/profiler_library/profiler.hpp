// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef COMP_LIBRARY_PROFILER_LIBRARY_PROFILER_HPP
#define COMP_LIBRARY_PROFILER_LIBRARY_PROFILER_HPP

#include "comp_library/profiler_library/events/time.hpp"

#ifdef DCA_NO_PROFILING
#include "comp_library/profiler_library/profilers/null_profiler.hpp"
#endif  // DCA_NO_PROFILING

#ifdef DCA_COUNTING_PROFILING
#include "comp_library/profiler_library/events/time_events.h"
#include "comp_library/profiler_library/profilers/counting_profiler.hpp"
#endif  // DCA_COUNTING_PROFILING

#ifdef DCA_PAPI_PROFILING
#include "comp_library/profiler_library/events/papi_events.h"
#include "comp_library/profiler_library/profilers/counting_profiler.hpp"
#endif  // DCA_PAPI_PROFILING

#ifdef DCA_CRAYPAT_PROFILING
#include "CraypatProfiler.h"
#endif  // DCA_CRAYPAT_PROFILING

#endif  // COMP_LIBRARY_PROFILER_LIBRARY_PROFILER_HPP
