// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef DCA_PROFILING_NULL_PROFILER_HPP
#define DCA_PROFILING_NULL_PROFILER_HPP

#include <string>

namespace dca {
namespace profiling {
// dca::profiling::

class NullProfiler {
public:
  NullProfiler(const char* /*functionName_*/, const char* /*fileName_*/, int /*line*/,
               bool /*bogusArgument*/ = true) {}
  NullProfiler(std::ostringstream& /*functionNameStrm*/, const char* /*fileName_*/, int /*line*/,
               bool /*bogusArgument*/ = true) {}

  static void start() {}
  static void stop(std::string /*fileName*/) {}
  template <typename ParallelProcessingType>
  static void stop(const ParallelProcessingType& /*parallelProcessing*/, std::string /*fileName*/) {}

  static void start_pthreading(int /*id*/) {}
  static void stop_pthreading(int /*id*/) {}
};

}  // profiling
}  // dca

#endif  // DCA_PROFILING_NULL_PROFILER_HPP
