// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: John Biddiscombe (john.biddiscombe@cscs.ch)
//
// A std::thread MC integrator that implements a threaded MC integration independent of the MC
// method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_CLUSTER_SOLVER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_CLUSTER_SOLVER_HPP

#include <atomic>
#include <iostream>
#include <future>
#include <stack>
#include <stdexcept>
#include <vector>

#include "dca/linalg/util/handle_functions.hpp"
#include "dca/parallel/stdthread/thread_pool/thread_pool.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/thread_task_handler.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/util/print_time.hpp"
#include "dca/parallel/util/get_workload.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <class QmciSolver>
class StdThreadQmciClusterSolver : public QmciSolver {
  using BaseClass = QmciSolver;
  using ThisType = StdThreadQmciClusterSolver<BaseClass>;

  using Data = typename BaseClass::DataType;
  using Parameters = typename BaseClass::ParametersType;
  using typename BaseClass::Concurrency;
  using typename BaseClass::Profiler;
  using typename BaseClass::Rng;

  using typename BaseClass::Walker;
  using typename BaseClass::Accumulator;
  using StdThreadAccumulatorType = stdthreadqmci::StdThreadQmciAccumulator<Accumulator>;

public:
  StdThreadQmciClusterSolver(Parameters& parameters_ref, Data& data_ref);

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

private:
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

  std::atomic<int> walk_finished_;
  std::atomic<int> measurements_remaining_;

  const int nr_walkers_;
  const int nr_accumulators_;

  ThreadTaskHandler thread_task_handler_;

  std::vector<Rng> rng_vector_;

  std::stack<StdThreadAccumulatorType*, std::vector<StdThreadAccumulatorType*>> accumulators_queue_;

  std::mutex mutex_merge_;
  std::mutex mutex_queue_;
  std::condition_variable queue_insertion_;
};

template <class QmciSolver>
StdThreadQmciClusterSolver<QmciSolver>::StdThreadQmciClusterSolver(Parameters& parameters_ref,
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

  // Create a sufficient amount of cublas handles, cuda streams and threads.
  linalg::util::resizeHandleContainer(thread_task_handler_.size());
  parallel::ThreadPool::get_instance().enlarge(thread_task_handler_.size());
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::initialize(int dca_iteration) {
  Profiler profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);

  BaseClass::initialize(dca_iteration);

  walk_finished_ = 0;
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::integrate() {
  Profiler profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);

  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "Threaded QMC integration has started: " << dca::util::print_time() << "\n"
              << std::endl;
  }

  measurements_remaining_ = parallel::util::getWorkload(parameters_.get_measurements(), concurrency_);

  if (concurrency_.id() == concurrency_.first())
    thread_task_handler_.print();

  std::vector<std::future<void>> futures;

  dca::profiling::WallTime start_time;

  auto& pool = dca::parallel::ThreadPool::get_instance();
  for (int i = 0; i < thread_task_handler_.size(); ++i) {
    if (thread_task_handler_.getTask(i) == "walker")
      futures.emplace_back(pool.enqueue(&ThisType::startWalker, this, i));
    else if (thread_task_handler_.getTask(i) == "accumulator")
      futures.emplace_back(pool.enqueue(&ThisType::startAccumulator, this, i));
    else if (thread_task_handler_.getTask(i) == "walker and accumulator")
      futures.emplace_back(pool.enqueue(&ThisType::startWalkerAndAccumulator, this, i));
    else
      throw std::logic_error("Thread task is undefined.");
  }

  for (auto& future : futures)
    future.get();
  assert(walk_finished_ == parameters_.get_walkers());

  dca::profiling::WallTime end_time;

  dca::profiling::Duration duration(end_time, start_time);
  total_time_ = duration.sec + 1.e-6 * duration.usec;

  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "Threaded on-node integration has ended: " << dca::util::print_time()
              << "\n\nTotal number of measurements: " << parameters_.get_measurements()
              << "\nQMC-time\t" << total_time_ << std::endl;
  }

  QmciSolver::accumulator_.finalize();
}

template <class QmciSolver>
template <typename dca_info_struct_t>
double StdThreadQmciClusterSolver<QmciSolver>::finalize(dca_info_struct_t& dca_info_struct) {
  Profiler profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);
  if (dca_iteration_ == parameters_.get_dca_iterations() - 1)
    BaseClass::computeErrorBars();

  double L2_Sigma_difference = QmciSolver::finalize(dca_info_struct);
  return L2_Sigma_difference;
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::startWalker(int id) {
  Profiler::start_threading(id);
  if (id == 0) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n\t\t QMCI starts\n" << std::endl;
  }

  const int walker_index = thread_task_handler_.walkerIDToRngIndex(id);
  Walker walker(parameters_, data_, rng_vector_[walker_index], id);

  walker.initialize();

  {
    Profiler profiler("thermalization", "stdthread-MC-walker", __LINE__, id);
    warmUp(walker, id);
  }

  StdThreadAccumulatorType* acc_ptr = nullptr;

  const int meas_to_do = parallel::util::getWorkload(
      parameters_.get_measurements(), parameters_.get_walkers(), walker_index, concurrency_);

  for (int meas_id = 0; meas_id < meas_to_do; ++meas_id) {
    {
      Profiler profiler("stdthread-MC-walker updating", "stdthread-MC-walker", __LINE__, id);
      walker.doSweep();
    }
    if (id == 0)
      walker.updateShell(meas_id, meas_to_do);

    {
      Profiler profiler("stdthread-MC-walker waiting", "stdthread-MC-walker", __LINE__, id);
      acc_ptr = nullptr;

      // Wait for available accumulators.
      {
        std::unique_lock<std::mutex> lock(mutex_queue_);
        queue_insertion_.wait(lock, [&]() { return !accumulators_queue_.empty(); });
        acc_ptr = accumulators_queue_.top();
        accumulators_queue_.pop();
      }

      acc_ptr->updateFrom(walker);
    }
  }

  // If this is the last walker signal to all the accumulators to exit the loop.
  if (++walk_finished_ == parameters_.get_walkers()) {
    std::lock_guard<std::mutex> lock(mutex_queue_);
    while (accumulators_queue_.size()) {
      accumulators_queue_.top()->notifyDone();
      accumulators_queue_.pop();
    }
  }

  if (id == 0 && concurrency_.id() == concurrency_.first()) {
    std::cout << "\n\t\t QMCI ends\n" << std::endl;
    walker.printSummary();
  }

  Profiler::stop_threading(id);
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::warmUp(Walker& walker, int id) {
  if (id == 0) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n\t\t warm-up starts\n" << std::endl;
  }

  for (int i = 0; i < parameters_.get_warm_up_sweeps(); i++) {
    walker.doSweep();

    if (id == 0)
      walker.updateShell(i, parameters_.get_warm_up_sweeps());
  }

  walker.is_thermalized() = true;

  if (id == 0) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n\t\t warm-up ends\n" << std::endl;
  }
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::startAccumulator(int id) {
  Profiler::start_threading(id);

  StdThreadAccumulatorType accumulator_obj(parameters_, data_, id);

  accumulator_obj.initialize(dca_iteration_);

  while (walk_finished_ < parameters_.get_walkers()) {
    {
      std::lock_guard<std::mutex> lock(mutex_queue_);
      accumulators_queue_.push(&accumulator_obj);
    }
    queue_insertion_.notify_one();

    {
      Profiler profiler("waiting", "stdthread-MC-accumulator", __LINE__, id);
      accumulator_obj.waitForQmciWalker();
    }

    {
      Profiler profiler("accumulating", "stdthread-MC-accumulator", __LINE__, id);
      accumulator_obj.measure();
    }
  }

  {
    std::lock_guard<std::mutex> lock(mutex_merge_);
    accumulator_obj.sumTo(QmciSolver::accumulator_);
  }

  Profiler::stop_threading(id);
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::startWalkerAndAccumulator(int id) {
  Profiler::start_threading(id);

  // Create and warm a walker.
  Walker walker(parameters_, data_, rng_vector_[id], id);
  walker.initialize();
  {
    Profiler profiler("thermalization", "stdthread-MC", __LINE__, id);
    warmUp(walker, id);
  }

  Accumulator accumulator_obj(parameters_, data_, id);
  accumulator_obj.initialize(dca_iteration_);

  const int n_meas = parallel::util::getWorkload(parameters_.get_measurements(),
                                                 parameters_.get_walkers(), id, concurrency_);

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

  ++walk_finished_;

  {
    std::lock_guard<std::mutex> lock(mutex_merge_);
    accumulator_obj.sumTo(QmciSolver::accumulator_);
  }

  Profiler::stop_threading(id);
}

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_CLUSTER_SOLVER_HPP
