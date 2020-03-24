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

#ifndef DCA_PARALLEL_UTIL_CALL_ONCE_PER_LOOP_HPP
#define DCA_PARALLEL_UTIL_CALL_ONCE_PER_LOOP_HPP

#include "dca/config/haves_defines.hpp"
#include "dca/config/threading.hpp"

namespace dca {
namespace util {
// dca::util::

struct OncePerLoopFlag {
  OncePerLoopFlag() : loop_done(-1) {}

  std::atomic<int> loop_done;
  dca::parallel::thread_traits::mutex_type mutex;
};

// This routine ensures f(args...) is called by a single thread for each value of loop index. Other
// threads calling this function with the same flag object and value of loop index wait on the
// completion of f(args...).
// Precondition: each call must use a non decreasing value of the loop index.
template <class F, class... Args>
void callOncePerLoop(OncePerLoopFlag& flag, const int loop_id, F&& f, Args&&... args)
{
    const int currently_done = flag.loop_done;

    if (loop_id < 0)
        throw(std::out_of_range("Negative loop index."));

    if (loop_id <= currently_done)
        return;
    else if (loop_id > currently_done + 1 && currently_done != -1)
        throw(std::logic_error("Loop id called out of order."));

    std::unique_lock<dca::parallel::thread_traits::mutex_type> lock(flag.mutex);
    // Check if flag.loop_done changed before locking the mutex.
    if (loop_id <= flag.loop_done)
        return;

    // Run the task.
    f(args...);

    flag.loop_done = loop_id;
}

}  // util
}  // dca

#endif  // DCA_PARALLEL_UTIL_CALL_ONCE_PER_LOOP_HPP
