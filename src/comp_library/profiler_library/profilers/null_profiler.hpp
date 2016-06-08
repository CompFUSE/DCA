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

#ifndef COMP_LIBRARY_PROFILER_LIBRARY_PROFILERS_NULL_PROFILER_HPP
#define COMP_LIBRARY_PROFILER_LIBRARY_PROFILERS_NULL_PROFILER_HPP

#include <string>

namespace PROFILER {

class NullProfiler {
public:
  NullProfiler(const char* /*functionName_*/, const char* /*fileName_*/, int /*line*/,
               bool /*bogusArgument*/ = true) {}

  NullProfiler(std::ostringstream& /*functionNameStrm*/, const char* /*fileName_*/, int /*line*/,
               bool /*bogusArgument*/ = true) {}

  static void start_pthreading(int /*id*/) {}

  static void stop_pthreading(int /*id*/) {}

  static void start() {}

  static void stop(std::string /*fileName*/) {}

  template <typename ParallelProcessingType>
  static void stop(const ParallelProcessingType& /*parallelProcessing*/, std::string /*fileName*/) {}
};
}

#endif  // COMP_LIBRARY_PROFILER_LIBRARY_PROFILERS_NULL_PROFILER_HPP
