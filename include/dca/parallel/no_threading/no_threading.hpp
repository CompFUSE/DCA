// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class is the equivalent of dca::parallel::stdthread for serial execution.

#ifndef DCA_PARALLEL_NO_THREADING_HPP
#define DCA_PARALLEL_NO_THREADING_HPP

#include <iostream>
#include <stdexcept>

namespace dca {
namespace parallel {
// dca::parallel::

class NoThreading {
public:
  // Execute the function f(id, num_tasks, args...) as num_tasks synchronous tasks with id in
  // [0, num_tasks - 1].
  template <class F, class... Args>
  void execute(int num_tasks, F&& f, Args&&... args) {
    for (int id = 0; id < num_tasks; ++id)
      f(id, num_tasks, args...);
  }

  // Returns \sum_{id = 0}^{num_threads -1} f(id, num_threads, args...).
  // Precondition: the return type of f can be initialized with 0.
  template <class F, class... Args>
  auto sumReduction(int num_threads, F&& f, Args&&... args) {
    using ReturnType = typename std::result_of<F(int, int, Args...)>::type;
    ReturnType result = 0;

    for (int id = 0; id < num_threads; ++id)
      result += f(id, 1, args...);

    return result;
  }

  friend std::ostream& operator<<(std::ostream& some_ostream, const NoThreading& this_concurrency);

private:
  constexpr static char parallel_type_str_[] = "NoThreading";
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_NO_THREADING_HPP
