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
#include <stdexcept>

#ifndef DCA_UTIL_CALL_ONCE_PER_LOOP_HPP
#define DCA_UTIL_CALL_ONCE_PER_LOOP_HPP

namespace dca {
namespace util {
// dca::util::

struct OncePerLoopFlag {
  OncePerLoopFlag() : loop_done(-1) {}

  bool increment_if_less_than(int val) {
    int old_loop = loop_done.load();
    do {
      if (old_loop >= val) {
        return false;
      }
      else if (old_loop < (val-1)) {
        throw std::logic_error("Loop id called out of order.");
      }
      else if (val < 0)
        throw(std::out_of_range("Negative loop index."));
    } while(!loop_done.compare_exchange_weak(old_loop, val));
    return true;
  }

  std::atomic<int> loop_done;
};

template <class F, class... Args>
void callOncePerLoop(OncePerLoopFlag& flag, const int loop_id, F&& f, Args&&... args) {

  if (!flag.increment_if_less_than(loop_id)) {
    return;
  }
  // Run the task.
  f(args...);
}

}  // util
}  // dca

#endif  // DCA_UTIL_CALL_ONCE_PER_LOOP_HPP
