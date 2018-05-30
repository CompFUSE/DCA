// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: John Biddiscombe (john.biddiscombe@cscs.ch)
//
// A std::thread jacket that implements a MC accumulator independent of the MC method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>

namespace dca {
namespace phys {
namespace solver {
namespace stdthreadqmci {
// dca::phys::solver::stdthreadqmci::

template <class QmciAccumulatorType>
class StdThreadQmciAccumulator : public QmciAccumulatorType {
  typedef StdThreadQmciAccumulator<QmciAccumulatorType> this_type;

public:
  template <class Parameters, class Data>
  StdThreadQmciAccumulator(Parameters& parameters_ref, Data& data_ref, int meas_to_do, int id);

  template <typename Walker>
  void updateFrom(Walker& walker);

  void waitForQmciWalker();

  void measure();

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sumTo(QmciAccumulatorType& other);

private:
  int thread_id_;
  int measurements_done_;
  int measurements_to_do_;
  bool measuring_;
  std::condition_variable start_measuring_;
  std::mutex mutex_accumulator_;
};

template <class QmciAccumulatorType>
template <class Parameters, class Data>
StdThreadQmciAccumulator<QmciAccumulatorType>::StdThreadQmciAccumulator(Parameters& parameters_ref,
                                                                        Data& data_ref,
                                                                        const int meas_to_do,
                                                                        const int id)
    : QmciAccumulatorType(parameters_ref, data_ref, id),
      thread_id_(id),
      measurements_done_(0),
      measurements_to_do_(meas_to_do),
      measuring_(false) {}

template <class QmciAccumulatorType>
template <typename Walker>
void StdThreadQmciAccumulator<QmciAccumulatorType>::updateFrom(Walker& walker) {
  {
    // take a lock and keep it until it goes out of scope
    std::unique_lock<std::mutex> lock(mutex_accumulator_);
    if (measuring_)
      throw std::logic_error(__FUNCTION__);

    QmciAccumulatorType::updateFrom(walker);
    measuring_ = true;

    if (thread_id_ == 1)
      walker.updateShell(measurements_done_, measurements_to_do_);
  }

  start_measuring_.notify_one();
}

template <class QmciAccumulatorType>
void StdThreadQmciAccumulator<QmciAccumulatorType>::waitForQmciWalker() {
  std::unique_lock<std::mutex> lock(mutex_accumulator_);
  start_measuring_.wait(lock, [this]() { return measuring_ == true; });
}

template <class QmciAccumulatorType>
void StdThreadQmciAccumulator<QmciAccumulatorType>::measure() {
  std::unique_lock<std::mutex> lock(mutex_accumulator_);
  QmciAccumulatorType::measure();
  measuring_ = false;
  ++measurements_done_;
}

template <class QmciAccumulatorType>
void StdThreadQmciAccumulator<QmciAccumulatorType>::sumTo(QmciAccumulatorType& other) {
  std::unique_lock<std::mutex> lock(mutex_accumulator_);
  QmciAccumulatorType::sumTo(other);
}

}  // stdthreadqmci
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP
