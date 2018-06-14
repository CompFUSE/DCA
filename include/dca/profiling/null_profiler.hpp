// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides a Null profiler to turn off profiling.

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

  static void start_threading(int /*id*/) {}
  static void stop_threading(int /*id*/) {}
};

}  // profiling
}  // dca

#endif  // DCA_PROFILING_NULL_PROFILER_HPP
