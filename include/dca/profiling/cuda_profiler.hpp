// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Annotate CPU functions in NVPROF. Depends on PTHREAD.

#ifndef DCA_PROFILING_CUDA_PROFILER_HPP
#define DCA_PROFILING_CUDA_PROFILER_HPP

#include <string>
#include <pthread.h>
#include <nvToolsExtCuda.h>

namespace dca {
namespace profiling {
// dca::profiling::

class CudaProfiler {
public:
  inline CudaProfiler(const std::string& functionName_, const std::string& fileName_, int line);

  inline CudaProfiler(const std::string& functionName_, const std::string& fileName_, int line,
                      int thread_id);

  inline ~CudaProfiler();

  static void start() {
    active_ = true;
  }

  static void stop() {
    active_ = false;
  }

  static void start_threading(int /*id*/){}
  static void stop_threading(int /*id*/){}
  template<class Concurrency>
  static void stop(Concurrency& /*conc*/, const std::string& /*name*/){}
  static void stop(const std::string& /*name*/){}

private:
  inline static bool active_ = false;
};

CudaProfiler::CudaProfiler(const std::string& function_name, const std::string& category_name,
                           int /*line*/) {
  if (active_) {
    nvtxNameOsThread(pthread_self(), "Master");
    nvtxRangePush((function_name + " - " + category_name).c_str());
  }
}

CudaProfiler::CudaProfiler(const std::string& function_name, const std::string& category_name,
                           int /*line*/, int id) {
  if (active_) {
    nvtxNameOsThread(pthread_self(), ("Thread" + std::to_string(id)).c_str());
    nvtxRangePush((function_name + " - " + category_name).c_str());
  }
}

CudaProfiler::~CudaProfiler() {
  if (active_)
    nvtxRangePop();
}

}  // namespace profiling
}  // namespace dca

#endif  // DCA_PROFILING_CUDA_PROFILER_HPP
