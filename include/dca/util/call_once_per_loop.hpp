// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Provides an utility to call a function only once per loop id.

#include <atomic>
#include <mutex>
#include <stdexcept>

#ifndef DCA_UTIL_CALL_ONCE_PER_LOOP_HPP
#define DCA_UTIL_CALL_ONCE_PER_LOOP_HPP

namespace dca {
namespace util {
// dca::util::

struct OncePerLoopFlag {
  OncePerLoopFlag() : loop_done(-1) {}

  std::atomic<int> loop_done;
  std::mutex mutex;
};

template <class F, class... Args>
void callOncePerLoop(OncePerLoopFlag& flag, const int loop_id, F&& f, Args&&... args) {
  if (loop_id < 0)
    throw(std::out_of_range("Negative loop index."));

  if (loop_id <= flag.loop_done)
    return;
  else if (loop_id > flag.loop_done + 1)
    throw(std::logic_error("Loop id called out of order."));

  std::unique_lock<std::mutex> lock(flag.mutex);
  if (loop_id == flag.loop_done)
    return;

  f(args...);

  ++flag.loop_done;
}

}  // util
}  // dca

#endif  // DCA_UTIL_CALL_ONCE_PER_LOOP_HPP
