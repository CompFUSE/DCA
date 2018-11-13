// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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

template <class QmciAccumulator>
class StdThreadQmciAccumulator : protected QmciAccumulator {
  using ThisType = StdThreadQmciAccumulator<QmciAccumulator>;
  using Parameters = typename QmciAccumulator::ParametersType;
  using Data = typename QmciAccumulator::DataType;

public:
  StdThreadQmciAccumulator(Parameters& parameters_ref, Data& data_ref, int meas_to_do, int id);

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

protected:
  using QmciAccumulator::get_Gflop;
  using QmciAccumulator::get_number_of_measurements;
  using QmciAccumulator::get_accumulated_sign;

private:
  using QmciAccumulator::data_;
  using QmciAccumulator::parameters_;

  int thread_id_;
  int measurements_done_;
  int measurements_to_do_;
  bool measuring_;
  std::condition_variable start_measuring_;
  std::mutex mutex_accumulator_;
};

template <class QmciAccumulator>
StdThreadQmciAccumulator<QmciAccumulator>::StdThreadQmciAccumulator(Parameters& parameters_ref,
                                                                    Data& data_ref,
                                                                    const int meas_to_do,
                                                                    const int id)
    : QmciAccumulator(parameters_ref, data_ref, id),
      thread_id_(id),
      measurements_done_(0),
      measurements_to_do_(meas_to_do),
      measuring_(false) {}

template <class QmciAccumulator>
StdThreadQmciAccumulator<QmciAccumulator>::~StdThreadQmciAccumulator() {}

template <class QmciAccumulator>
template <typename walker_type>
void StdThreadQmciAccumulator<QmciAccumulator>::updateFrom(walker_type& walker) {
  {
    // take a lock and keep it until it goes out of scope
    std::unique_lock<std::mutex> lock(mutex_accumulator_);
    if (measuring_)
      throw std::logic_error(__FUNCTION__);

    QmciAccumulator::updateFrom(walker);
    measuring_ = true;

    if (thread_id_ == 1)
      walker.updateShell(measurements_done_, measurements_to_do_);
  }

  start_measuring_.notify_one();
}

template <class QmciAccumulator>
void StdThreadQmciAccumulator<QmciAccumulator>::waitForQmciWalker() {
  std::unique_lock<std::mutex> lock(mutex_accumulator_);
  start_measuring_.wait(lock, [this]() { return measuring_ == true; });
}

template <class QmciAccumulator>
void StdThreadQmciAccumulator<QmciAccumulator>::measure() {
  std::unique_lock<std::mutex> lock(mutex_accumulator_);
  QmciAccumulator::measure();
  measuring_ = false;
  ++measurements_done_;
}

template <class QmciAccumulator>
void StdThreadQmciAccumulator<QmciAccumulator>::sumTo(QmciAccumulator& other) {
  std::unique_lock<std::mutex> lock(mutex_accumulator_);
  QmciAccumulator::sumTo(other);
}

}  // stdthreadqmci
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP
