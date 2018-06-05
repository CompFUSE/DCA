// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: John Biddiscombe (john.biddiscombe@cscs.ch)
//
// A std::thread MC integrator that implements a threaded MC integration independent of the MC
// method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_CLUSTER_SOLVER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_CLUSTER_SOLVER_HPP

#include <atomic>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/thread_task_handler.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/util/print_time.hpp"
#include "dca/parallel/util/get_workload.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <class QmcSolver>
class StdThreadQmciClusterSolver : public QmcSolver {
  using BaseClass = QmcSolver;
  using ThisType = StdThreadQmciClusterSolver<BaseClass>;
  using Data = typename BaseClass::DataType;
  using Parameters = typename BaseClass::ParametersType;
  using typename BaseClass::Concurrency;
  using typename BaseClass::Profiler;
  using typename BaseClass::Rng;

  using typename BaseClass::Walker;
  using typename BaseClass::Accumulator;
  using StdThreadAccumulatorType = stdthreadqmci::StdThreadQmciAccumulator<Accumulator>;

  using BaseClass::instantiateWalker;

  using Pair = std::pair<ThisType*, int>;

public:
  StdThreadQmciClusterSolver(Parameters& parameters_ref, Data& data_ref);

  template <typename Writer>
  void write(Writer& reader);

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

private:
  static void* startWalkerStatic(Pair data);
  static void* startAccumulatorStatic(Pair data);
  static void* startWalkerAndAccumulatorStatic(Pair data);

  void startWalker(int id);
  void startAccumulator(int id);
  void startWalkerAndAccumulator(int id);

  void warmUp(Walker& walker, int id);

private:
  using BaseClass::parameters_;
  using BaseClass::data_;
  using BaseClass::concurrency_;
  using BaseClass::total_time_;
  using BaseClass::dca_iteration_;
  using BaseClass::accumulator_;

  std::atomic<int> acc_finished_;

  std::atomic<int> acc_finished;
  std::atomic<int> measurements_remaining_;

  const int nr_walkers_;
  const int nr_accumulators_;

  ThreadTaskHandler thread_task_handler_;

  std::vector<Rng> rng_vector_;

  std::queue<StdThreadAccumulatorType*> accumulators_queue_;

  std::mutex mutex_merge_;
  std::mutex mutex_queue_;
};

template <class QmcSolver>
StdThreadQmciClusterSolver<QmcSolver>::StdThreadQmciClusterSolver(Parameters& parameters_ref,
                                                                  Data& data_ref)
    : BaseClass(parameters_ref, data_ref),

      nr_walkers_(parameters_.get_walkers()),
      nr_accumulators_(parameters_.get_accumulators()),

      thread_task_handler_(nr_walkers_, nr_accumulators_,
                           parameters_ref.shared_walk_and_accumulation_thread()),

      accumulators_queue_() {
  if (nr_walkers_ < 1 || nr_accumulators_ < 1) {
    throw std::logic_error(
        "Both the number of walkers and the number of accumulators must be at least 1.");
  }

  for (int i = 0; i < nr_walkers_; ++i) {
    rng_vector_.emplace_back(concurrency_.id(), concurrency_.number_of_processors(),
                             parameters_.get_seed());
  }
}

template <class QmcSolver>
template <typename Writer>
void StdThreadQmciClusterSolver<QmcSolver>::write(Writer& writer) {
  BaseClass::write(writer);
  // accumulator.write(writer);
}

template <class QmcSolver>
void StdThreadQmciClusterSolver<QmcSolver>::initialize(int dca_iteration) {
  Profiler profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);

  BaseClass::initialize(dca_iteration);

  acc_finished_ = 0;
}

template <class QmcSolver>
void StdThreadQmciClusterSolver<QmcSolver>::integrate() {
  Profiler profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);

  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "Threaded QMC integration has started: " << dca::util::print_time() << "\n"
              << std::endl;
  }

  measurements_remaining_ =
      parallel::util::getWorkload(parameters_.get_measurements(), 1, 0, concurrency_);

  std::vector<std::thread> threads;
  std::vector<Pair> data;

  {
    if (concurrency_.id() == concurrency_.first())
      thread_task_handler_.print();

    dca::profiling::WallTime start_time;

    for (int i = 0; i < thread_task_handler_.size(); ++i) {
      data.push_back(std::make_pair(this, i));

      if (thread_task_handler_.getTask(i) == "walker")
        threads.push_back(std::thread(startWalkerStatic, data.back()));
      else if (thread_task_handler_.getTask(i) == "accumulator")
        threads.push_back(std::thread(startAccumulatorStatic, data.back()));
      else if (thread_task_handler_.getTask(i) == "walker and accumulator")
        threads.push_back(std::thread(startWalkerAndAccumulatorStatic, data.back()));
      else
        throw std::logic_error("Thread task is undefined.");
    }

    for (auto& thread : threads)
      thread.join();

    dca::profiling::WallTime end_time;

    dca::profiling::Duration duration(end_time, start_time);
    total_time_ = duration.sec + 1.e-6 * duration.usec;
  }

  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "Threaded on-node integration has ended: " << dca::util::print_time()
              << "\n\nTotal number of measurements: " << parameters_.get_measurements()
              << "\nQMC-time\t" << total_time_ << std::endl;
  }

  accumulator_.finalize();
}

template <class QmcSolver>
template <typename dca_info_struct_t>
double StdThreadQmciClusterSolver<QmcSolver>::finalize(dca_info_struct_t& dca_info_struct) {
  Profiler profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);
  if (dca_iteration_ == parameters_.get_dca_iterations() - 1)
    BaseClass::computeErrorBars();

  double L2_Sigma_difference = BaseClass::finalize(dca_info_struct);
  return L2_Sigma_difference;
}

template <class QmcSolver>
void* StdThreadQmciClusterSolver<QmcSolver>::startWalkerStatic(Pair data) {
  Profiler::start_threading(data.second);

  data.first->startWalker(data.second);

  Profiler::stop_threading(data.second);

  return NULL;
}

template <class QmcSolver>
void* StdThreadQmciClusterSolver<QmcSolver>::startAccumulatorStatic(Pair data) {
  Profiler::start_threading(data.second);

  data.first->startAccumulator(data.second);

  Profiler::stop_threading(data.second);

  return NULL;
}

template <class QmcSolver>
void* StdThreadQmciClusterSolver<QmcSolver>::startWalkerAndAccumulatorStatic(Pair data) {
  Profiler::start_threading(data.second);

  data.first->startWalkerAndAccumulator(data.second);

  Profiler::stop_threading(data.second);

  return NULL;
}

template <class QmcSolver>
void StdThreadQmciClusterSolver<QmcSolver>::startWalker(int id) {
  if (id == 0) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n\t\t QMCI starts\n" << std::endl;
  }

  const int rng_index = thread_task_handler_.walkerIDToRngIndex(id);
  Walker walker(instantiateWalker(rng_vector_[rng_index], id));

  walker.initialize();

  {
    Profiler profiler("thermalization", "stdthread-MC-walker", __LINE__, id);
    warmUp(walker, id);
  }

  StdThreadAccumulatorType* acc_ptr(nullptr);

  while (--measurements_remaining_ >= 0) {
    {
      Profiler profiler("stdthread-MC-walker updating", "stdthread-MC-walker", __LINE__, id);
      walker.doSweep();
    }

    {
      Profiler profiler("stdthread-MC-walker waiting", "stdthread-MC-walker", __LINE__, id);
      acc_ptr = nullptr;

      while (acc_ptr == nullptr) {  // checking for available accumulators
        std::unique_lock<std::mutex> lock(mutex_queue_);
        if (!accumulators_queue_.empty()) {
          acc_ptr = accumulators_queue_.front();
          accumulators_queue_.pop();
        }
      }

      acc_ptr->updateFrom(walker);
    }
  }

  //  #ifdef DCA_WITH_QMC_BIT
  //  mutex_numerical_error.lock();
  //  accumulator.get_error_distribution() += walker.get_error_distribution();
  //  mutex_numerical_error.unlock();
  //  #endif  // DCA_WITH_QMC_BIT

  if (id == 0) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n\t\t QMCI ends\n" << std::endl;
  }
}

template <class QmcSolver>
void StdThreadQmciClusterSolver<QmcSolver>::warmUp(Walker& walker, int id) {
  if (id == 0) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n\t\t warm-up starts\n" << std::endl;
  }

  for (int i = 0; i < parameters_.get_warm_up_sweeps(); i++) {
    walker.doSweep();

    if (id == 0)
      walker.updateShell(i, parameters_.get_warm_up_sweeps());
  }

  walker.markThermalized();

  if (id == 0) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n\t\t warm-up ends\n" << std::endl;
  }
}

template <class qmci_integrator_type>
void StdThreadQmciClusterSolver<qmci_integrator_type>::startAccumulator(int id) {
  const int n_meas =
      parallel::util::getWorkload(parameters_.get_measurements(), parameters_.get_accumulators(),
                                  thread_task_handler_.IDToAccumIndex(id), concurrency_);

  StdThreadAccumulatorType accumulator_obj(parameters_, data_, n_meas, id);

  accumulator_obj.initialize(dca_iteration_);

  for (int i = 0; i < n_meas; ++i) {
    {
      std::lock_guard<std::mutex> lock(mutex_queue_);
      accumulators_queue_.push(&accumulator_obj);
    }

    {
      Profiler profiler("stdthread-accumulator waiting", "stdthread-MC-accumulator", __LINE__, id);
      accumulator_obj.waitForQmciWalker();
    }

    {
      Profiler profiler("stdthread-accumulator accumulating", "stdthread-MC-accumulator", __LINE__,
                        id);
      accumulator_obj.measure();
    }
  }

  ++acc_finished_;
  {
    std::lock_guard<std::mutex> lock(mutex_merge_);
    accumulator_obj.sumTo(accumulator_);
  }
}

template <class QmcSolver>
void StdThreadQmciClusterSolver<QmcSolver>::startWalkerAndAccumulator(int id) {
  // Create and warm a walker.
  Walker walker(instantiateWalker(rng_vector_[id], id));
  walker.initialize();
  {
    Profiler profiler("thermalization", "stdthread-MC", __LINE__, id);
    warmUp(walker, id);
  }

  Accumulator accumulator_obj(parameters_, data_, id);
  accumulator_obj.initialize(dca_iteration_);

  const int n_meas = parallel::util::getWorkload(parameters_.get_measurements(),
                                                 parameters_.get_accumulators(), id, concurrency_);

  for (int i = 0; i < n_meas; ++i) {
    {
      Profiler profiler("Walker updating", "stdthread-MC", __LINE__, id);
      walker.doSweep();
    }
    {
      Profiler profiler("Accumulator measuring", "stdthread-MC", __LINE__, id);
      accumulator_obj.updateFrom(walker);
      accumulator_obj.measure();
    }
    if (id == 0)
      walker.updateShell(i, n_meas);
  }

  ++acc_finished_;
  std::lock_guard<std::mutex> lock(mutex_merge_);
  accumulator_obj.sumTo(accumulator_);
}

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_CLUSTER_SOLVER_HPP
