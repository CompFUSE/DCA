// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// A std::thread jacket that implements a MC accumulator independent of the MC method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <stdexcept>
#include <atomic>

namespace dca {
namespace phys {
namespace solver {
namespace stdthreadqmci {

template <class qmci_accumulator_type>
class stdthread_qmci_accumulator : protected qmci_accumulator_type {
  typedef typename qmci_accumulator_type::my_parameters_type parameters_type;
  typedef typename qmci_accumulator_type::my_MOMS_type MOMS_type;

  typedef stdthread_qmci_accumulator<qmci_accumulator_type> this_type;

public:
  stdthread_qmci_accumulator(parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id);

  ~stdthread_qmci_accumulator();

  using qmci_accumulator_type::finalize;
  using qmci_accumulator_type::initialize;
  // using qmci_accumulator_type::to_JSON;
  using qmci_accumulator_type::get_configuration;

  template <typename walker_type>
  void update_from(walker_type& walker);

  void wait_for_qmci_walker();

  void measure(std::mutex& mutex_queue, std::queue<this_type*>& accumulators_queue);

  // void sum_to(qmci_accumulator_type& accumulator_obj);
  // int get_expansion_order();

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sum_to(qmci_accumulator_type& other);

protected:
  using qmci_accumulator_type::get_Gflop;
  using qmci_accumulator_type::get_number_of_measurements;
  using qmci_accumulator_type::get_sign;

private:
  using qmci_accumulator_type::MOMS;
  using qmci_accumulator_type::parameters;

  int thread_id;
  std::atomic<bool> measuring;
  std::condition_variable start_measuring;
  std::mutex mutex_accumulator;
};

template <class qmci_accumulator_type>
stdthread_qmci_accumulator<qmci_accumulator_type>::stdthread_qmci_accumulator(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id)
    : qmci_accumulator_type(parameters_ref, MOMS_ref, id), thread_id(id), measuring(false) {}

template <class qmci_accumulator_type>
stdthread_qmci_accumulator<qmci_accumulator_type>::~stdthread_qmci_accumulator() {}

template <class qmci_accumulator_type>
template <typename walker_type>
void stdthread_qmci_accumulator<qmci_accumulator_type>::update_from(walker_type& walker) {
  {
    // take a lock and keep it until it goes out of scope
    {
      std::lock_guard<std::mutex> lock(mutex_accumulator);
      if (measuring)
        throw std::logic_error(__FUNCTION__);

      qmci_accumulator_type::update_from(walker);
      measuring = true;
    }
    start_measuring.notify_one();
  }
}

template <class qmci_accumulator_type>
void stdthread_qmci_accumulator<qmci_accumulator_type>::wait_for_qmci_walker() {
  std::unique_lock<std::mutex> lock(mutex_accumulator);
  start_measuring.wait(lock, [this]() { return measuring==true; });
}

template <class qmci_accumulator_type>
void stdthread_qmci_accumulator<qmci_accumulator_type>::measure(
    std::mutex& /*mutex_queue*/, std::queue<this_type*>& /*accumulators_queue*/)
{
  std::unique_lock<std::mutex> lock(mutex_accumulator);
  qmci_accumulator_type::measure();
  measuring = false;
}

template <class qmci_accumulator_type>
void stdthread_qmci_accumulator<qmci_accumulator_type>::sum_to(qmci_accumulator_type& other) {
  std::unique_lock<std::mutex> lock(mutex_accumulator);
  qmci_accumulator_type::sum_to(other);
}

}  // namespace stdthreadqmci
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP
