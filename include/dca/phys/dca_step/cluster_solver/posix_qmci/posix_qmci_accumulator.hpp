// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// A posix jacket that implements a MC accumulator independent of the MC method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_POSIX_QMCI_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_POSIX_QMCI_ACCUMULATOR_HPP

#include <stdexcept>

#include "dca/parallel/pthreading/lock.hpp"
#include "dca/parallel/pthreading/pthreading.hpp"
#include "dca/parallel/pthreading/conditional_variable.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace posixqmci {
// dca::phys::solver::posixqmci::

template <class QmciAccumulatorType>
class PosixQmciAccumulator : public QmciAccumulatorType {
  using this_type = PosixQmciAccumulator<QmciAccumulatorType>;

public:
  PosixQmciAccumulator(QmciAccumulatorType&& accumulator);

  // Copies the state of a walker into the accumulator. It should be called by the walker thread.
  template <class Walker>
  void updateFrom(Walker& walker);

  // Blocks the accumulation thread until updateFrom is called by a walker thread.
  void waitForQmciWalker();

  // Accumulate the stored configuration. Should be called by an accumulator thread.
  void measure();

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sumTo(QmciAccumulatorType& other);

private:
  bool measuring_;
  parallel::ConditionalVariable start_measuring_cond_;
  parallel::Lock<parallel::Pthreading> accumulator_lock_;
};

template <class QmciAccumulatorType>
PosixQmciAccumulator<QmciAccumulatorType>::PosixQmciAccumulator(QmciAccumulatorType&& accumulator)
    : QmciAccumulatorType(std::move(accumulator)), measuring_(false) {}

template <class QmciAccumulatorType>
template <typename Walker>
void PosixQmciAccumulator<QmciAccumulatorType>::updateFrom(Walker& walker) {
  {
    accumulator_lock_.lock();

    if (measuring_)
      throw std::logic_error(__FUNCTION__);

    QmciAccumulatorType::updateFrom(walker);

    measuring_ = true;

    // Unblocks the accumulation thread.
    start_measuring_cond_.signal();

    accumulator_lock_.unlock();
  }
}

template <class QmciAccumulatorType>
void PosixQmciAccumulator<QmciAccumulatorType>::waitForQmciWalker() {
  accumulator_lock_.lock();

  while (!measuring_)
    start_measuring_cond_.wait(accumulator_lock_);

  accumulator_lock_.unlock();
}

template <class QmciAccumulatorType>
void PosixQmciAccumulator<QmciAccumulatorType>::measure() {
  accumulator_lock_.lock();

  QmciAccumulatorType::measure();

  measuring_ = false;

  accumulator_lock_.unlock();
}

template <class QmciAccumulatorType>
void PosixQmciAccumulator<QmciAccumulatorType>::sumTo(QmciAccumulatorType& other) {
  accumulator_lock_.lock();
  if (measuring_)
    throw std::logic_error(__FUNCTION__);
  QmciAccumulatorType::sumTo(other);
  accumulator_lock_.unlock();
}

}  // posixqmci
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_POSIX_QMCI_ACCUMULATOR_HPP
