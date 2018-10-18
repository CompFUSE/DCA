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

template <class qmci_accumulator_type>
class stdthread_qmci_accumulator : protected qmci_accumulator_type {
  using ThisType = stdthread_qmci_accumulator<qmci_accumulator_type>;
  using Parameters = typename qmci_accumulator_type::my_parameters_type;
  using Data = typename qmci_accumulator_type::DataType;

public:
  stdthread_qmci_accumulator(Parameters& parameters_ref, Data& data_ref, int meas_to_do, int id);

  ~stdthread_qmci_accumulator();

  using qmci_accumulator_type::finalize;
  using qmci_accumulator_type::initialize;
  // using qmci_accumulator_type::to_JSON;
  using qmci_accumulator_type::get_configuration;

  template <typename walker_type>
  void update_from(walker_type& walker);

  void wait_for_qmci_walker();

  void measure();

  // int get_expansion_order();

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sum_to(qmci_accumulator_type& other);

protected:
  using qmci_accumulator_type::get_Gflop;
  using qmci_accumulator_type::get_number_of_measurements;
  using qmci_accumulator_type::get_accumulated_sign;

private:
  using qmci_accumulator_type::data_;
  using qmci_accumulator_type::parameters;

  int thread_id;
  int measurements_done_;
  int measurements_to_do_;
  bool measuring;
  std::condition_variable start_measuring;
  std::mutex mutex_accumulator;
};

template <class qmci_accumulator_type>
stdthread_qmci_accumulator<qmci_accumulator_type>::stdthread_qmci_accumulator(
    Parameters& parameters_ref, Data& data_ref, const int meas_to_do, const int id)
    : qmci_accumulator_type(parameters_ref, data_ref, id),
      thread_id(id),
      measurements_done_(0),
      measurements_to_do_(meas_to_do),
      measuring(false) {}

template <class qmci_accumulator_type>
stdthread_qmci_accumulator<qmci_accumulator_type>::~stdthread_qmci_accumulator() {}

template <class qmci_accumulator_type>
template <typename walker_type>
void stdthread_qmci_accumulator<qmci_accumulator_type>::update_from(walker_type& walker) {
  {
    // take a lock and keep it until it goes out of scope
    std::unique_lock<std::mutex> lock(mutex_accumulator);
    if (measuring)
      throw std::logic_error(__FUNCTION__);

    qmci_accumulator_type::update_from(walker);
    measuring = true;

    if (thread_id == 1)
      walker.update_shell(measurements_done_, measurements_to_do_);
  }

  start_measuring.notify_one();
}

template <class qmci_accumulator_type>
void stdthread_qmci_accumulator<qmci_accumulator_type>::wait_for_qmci_walker() {
  std::unique_lock<std::mutex> lock(mutex_accumulator);
  start_measuring.wait(lock, [this]() { return measuring == true; });
}

template <class qmci_accumulator_type>
void stdthread_qmci_accumulator<qmci_accumulator_type>::measure() {
  std::unique_lock<std::mutex> lock(mutex_accumulator);
  qmci_accumulator_type::measure();
  measuring = false;
  ++measurements_done_;
}

template <class qmci_accumulator_type>
void stdthread_qmci_accumulator<qmci_accumulator_type>::sum_to(qmci_accumulator_type& other) {
  std::unique_lock<std::mutex> lock(mutex_accumulator);
  qmci_accumulator_type::sum_to(other);
}

}  // stdthreadqmci
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP
