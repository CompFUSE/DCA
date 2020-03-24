// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: John Biddiscombe (john.biddiscombe@cscs.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// A std::thread jacket that implements a MC accumulator independent of the MC method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP

#include <atomic>
#include <queue>
#include <stdexcept>

#include "dca/config/threading.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace stdthreadqmci {
// dca::phys::solver::stdthreadqmci::

template <class QmciAccumulator>
class StdThreadQmciAccumulator : public QmciAccumulator {
  using ThisType = StdThreadQmciAccumulator<QmciAccumulator>;
  using Parameters = typename QmciAccumulator::ParametersType;
  using Data = typename QmciAccumulator::DataType;

public:
  StdThreadQmciAccumulator(Parameters& parameters_ref, Data& data_ref, int id);

  ~StdThreadQmciAccumulator();

  using QmciAccumulator::finalize;
  using QmciAccumulator::initialize;
  using QmciAccumulator::get_configuration;

  template <typename walker_type>
  void updateFrom(walker_type& walker);

  void waitForQmciWalker();

  void measure();

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sumTo(QmciAccumulator& other);

  // Signals that this object will not need to perform any more accumulation.
  void notifyDone();

  bool done() const {
    return done_;
  }

protected:
  using QmciAccumulator::get_Gflop;
  using QmciAccumulator::get_number_of_measurements;
  using QmciAccumulator::get_accumulated_sign;

private:
  using QmciAccumulator::data_;
  using QmciAccumulator::parameters_;

  int thread_id_;
  bool measuring_;
  std::atomic<bool> done_;
  dca::parallel::thread_traits::condition_variable_type start_measuring_;
  dca::parallel::thread_traits::mutex_type mutex_accumulator_;
};

template <class QmciAccumulator>
StdThreadQmciAccumulator<QmciAccumulator>::StdThreadQmciAccumulator(Parameters& parameters_ref,
                                                                    Data& data_ref, const int id)
    : QmciAccumulator(parameters_ref, data_ref, id), thread_id_(id), measuring_(false), done_(false) {}

template <class QmciAccumulator>
StdThreadQmciAccumulator<QmciAccumulator>::~StdThreadQmciAccumulator() {}

template <class QmciAccumulator>
template <typename walker_type>
void StdThreadQmciAccumulator<QmciAccumulator>::updateFrom(walker_type& walker) {
  {
    // take a lock and keep it until it goes out of scope
    dca::parallel::thread_traits::unique_lock lock(mutex_accumulator_);
    if (measuring_)
      throw std::logic_error(__FUNCTION__);

    QmciAccumulator::updateFrom(walker);
    measuring_ = true;
  }

  start_measuring_.notify_one();
}

template <class QmciAccumulator>
void StdThreadQmciAccumulator<QmciAccumulator>::waitForQmciWalker() {
  dca::parallel::thread_traits::unique_lock lock(mutex_accumulator_);
  start_measuring_.wait(lock, [this]() { return measuring_ || done_; });
}

template <class QmciAccumulator>
void StdThreadQmciAccumulator<QmciAccumulator>::measure() {
  dca::parallel::thread_traits::scoped_lock lock(mutex_accumulator_);

  if (done_)
    return;
  assert(measuring_);

  QmciAccumulator::measure();
  measuring_ = false;
}

template <class QmciAccumulator>
void StdThreadQmciAccumulator<QmciAccumulator>::sumTo(QmciAccumulator& other) {
  dca::parallel::thread_traits::scoped_lock lock(mutex_accumulator_);
  QmciAccumulator::sumTo(other);
}

template <class QmciAccumulator>
void StdThreadQmciAccumulator<QmciAccumulator>::notifyDone() {
  done_ = true;
  start_measuring_.notify_one();
}

}  // stdthreadqmci
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP
