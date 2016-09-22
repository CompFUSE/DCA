// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Peter Staar (taa@zurich.ibm.com)
//
// TODO: This class seems quite redundant.

#ifndef DCA_PARALLEL_THREAD_MANAGER_SUM_HPP
#define DCA_PARALLEL_THREAD_MANAGER_SUM_HPP

#include <utility>

namespace dca {
namespace concurrency {
// dca::concurrency::

template <typename Concurrency>
class ThreadManagerSum {
public:
  ThreadManagerSum(const Concurrency& concurrency) : concurrency_(concurrency) {}

  // TODO: Add const to function parameter 'dmn'.
  template <typename Domain>
  std::pair<int, int> get_bounds(Domain& dmn) const {
    return concurrency_.get_bounds(dmn);
  }

  template <typename T>
  bool sum_and_check(T& result) const {
    concurrency_.sum(result);
    return true;
  }

private:
  const Concurrency& concurrency_;
};

}  // concurrency
}  // dca

#endif  // DCA_PARALLEL_THREAD_MANAGER_SUM_HPP
