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
#include <queue>
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

template <class qmci_integrator_type>
class StdThreadQmciClusterSolver : public qmci_integrator_type {
  using ThisType = StdThreadQmciClusterSolver<qmci_integrator_type>;

  using Data = typename qmci_integrator_type::DataType;
  typedef typename qmci_integrator_type::this_parameters_type parameters_type;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  using random_number_generator = typename parameters_type::random_number_generator;

  typedef typename qmci_integrator_type::walker_type walker_type;
  typedef typename qmci_integrator_type::accumulator_type accumulator_type;

  typedef StdThreadQmciClusterSolver<qmci_integrator_type> this_type;
  typedef stdthreadqmci::stdthread_qmci_accumulator<accumulator_type> stdthread_accumulator_type;

  typedef std::pair<this_type*, int> pair_type;

public:
  StdThreadQmciClusterSolver(parameters_type& parameters_ref, Data& data_ref);

  template <typename Writer>
  void write(Writer& reader);

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

private:
  void start_walker(int id);
  void start_accumulator(int id);
  void start_walker_and_accumulator(int id);

  void warm_up(walker_type& walker, int id);

  // TODO: Are the following using statements redundant and can therefore be removed?
  using qmci_integrator_type::compute_error_bars;
  using qmci_integrator_type::symmetrize_measurements;

private:
  using qmci_integrator_type::parameters;
  using qmci_integrator_type::data_;
  using qmci_integrator_type::concurrency;

  using qmci_integrator_type::total_time;

  using qmci_integrator_type::DCA_iteration;

  using qmci_integrator_type::accumulator;

  std::atomic<int> acc_finished;
  std::atomic<int> measurements_remaining_;

  const int nr_walkers;
  const int nr_accumulators;

  ThreadTaskHandler thread_task_handler_;

  std::vector<random_number_generator> rng_vector;

  std::queue<stdthread_accumulator_type*> accumulators_queue;

  std::mutex mutex_print;
  std::mutex mutex_merge;
  std::mutex mutex_queue;
  std::mutex mutex_acc_finished;
  std::mutex mutex_numerical_error;
  std::condition_variable queue_insertion_;
};

template <class qmci_integrator_type>
StdThreadQmciClusterSolver<qmci_integrator_type>::StdThreadQmciClusterSolver(
    parameters_type& parameters_ref, Data& data_ref)
    : qmci_integrator_type(parameters_ref, data_ref),

      nr_walkers(parameters.get_walkers()),
      nr_accumulators(parameters.get_accumulators()),

      thread_task_handler_(nr_walkers, nr_accumulators,
                           parameters.shared_walk_and_accumulation_thread()),

      accumulators_queue() {
  if (nr_walkers < 1 || nr_accumulators < 1) {
    throw std::logic_error(
        "Both the number of walkers and the number of accumulators must be at least 1.");
  }

  for (int i = 0; i < nr_walkers; ++i) {
    rng_vector.emplace_back(concurrency.id(), concurrency.number_of_processors(),
                            parameters.get_seed());
  }

  // Create a sufficient amount of cublas handles, cuda streams and threads.
  linalg::util::resizeHandleContainer(thread_task_handler_.size());
  parallel::ThreadPool::get_instance().enlarge(thread_task_handler_.size());
}

template <class qmci_integrator_type>
template <typename Writer>
void StdThreadQmciClusterSolver<qmci_integrator_type>::write(Writer& writer) {
  qmci_integrator_type::write(writer);
  // accumulator.write(writer);
}

template <class qmci_integrator_type>
void StdThreadQmciClusterSolver<qmci_integrator_type>::initialize(int dca_iteration) {
  profiler_type profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);

  qmci_integrator_type::initialize(dca_iteration);

  acc_finished = 0;
}

template <class qmci_integrator_type>
void StdThreadQmciClusterSolver<qmci_integrator_type>::integrate() {
  profiler_type profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);

  if (concurrency.id() == concurrency.first()) {
    std::cout << "Threaded QMC integration has started: " << dca::util::print_time() << "\n"
              << std::endl;
  }

  measurements_remaining_ = parallel::util::getWorkload(parameters.get_measurements(), concurrency);

  if (concurrency.id() == concurrency.first())
    thread_task_handler_.print();

  std::vector<std::future<void>> futures;

  dca::profiling::WallTime start_time;

  auto& pool = dca::parallel::ThreadPool::get_instance();
  for (int i = 0; i < thread_task_handler_.size(); ++i) {
    if (thread_task_handler_.getTask(i) == "walker")
      futures.emplace_back(pool.enqueue(&ThisType::start_walker, this, i));
    else if (thread_task_handler_.getTask(i) == "accumulator")
      futures.emplace_back(pool.enqueue(&ThisType::start_accumulator, this, i));
    else if (thread_task_handler_.getTask(i) == "walker and accumulator")
      futures.emplace_back(pool.enqueue(&ThisType::start_walker_and_accumulator, this, i));
    else
      throw std::logic_error("Thread task is undefined.");
  }

  for (auto& future : futures)
    future.get();

  dca::profiling::WallTime end_time;

  dca::profiling::Duration duration(end_time, start_time);
  total_time = duration.sec + 1.e-6 * duration.usec;

  if (concurrency.id() == concurrency.first()) {
    std::cout << "Threaded on-node integration has ended: " << dca::util::print_time()
              << "\n\nTotal number of measurements: " << parameters.get_measurements()
              << "\nQMC-time\t" << total_time << std::endl;
  }

  qmci_integrator_type::accumulator.finalize();
}

template <class qmci_integrator_type>
template <typename dca_info_struct_t>
double StdThreadQmciClusterSolver<qmci_integrator_type>::finalize(dca_info_struct_t& dca_info_struct) {
  profiler_type profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);
  if (DCA_iteration == parameters.get_dca_iterations() - 1)
    compute_error_bars();

  double L2_Sigma_difference = qmci_integrator_type::finalize(dca_info_struct);
  return L2_Sigma_difference;
}

template <class qmci_integrator_type>
void StdThreadQmciClusterSolver<qmci_integrator_type>::start_walker(int id) {
  profiler_type::start_threading(id);
  if (id == 0) {
    if (concurrency.id() == concurrency.first())
      std::cout << "\n\t\t QMCI starts\n" << std::endl;
  }

  const int rng_index = thread_task_handler_.walkerIDToRngIndex(id);
  walker_type walker(parameters, data_, rng_vector[rng_index], id);

  walker.initialize();

  {
    profiler_type profiler("thermalization", "stdthread-MC-walker", __LINE__, id);
    warm_up(walker, id);
  }

  stdthread_accumulator_type* acc_ptr(NULL);

  while (--measurements_remaining_ >= 0) {
    {
      profiler_type profiler("stdthread-MC-walker updating", "stdthread-MC-walker", __LINE__, id);
      walker.do_sweep();
    }

    {
      profiler_type profiler("stdthread-MC-walker waiting", "stdthread-MC-walker", __LINE__, id);
      acc_ptr = nullptr;

      // Wait for available accumulators.
      {
        std::unique_lock<std::mutex> lock(mutex_queue);
        queue_insertion_.wait(lock, [&]() { return !accumulators_queue.empty(); });
        acc_ptr = accumulators_queue.front();
        accumulators_queue.pop();
      }

      acc_ptr->update_from(walker);
    }
  }

  //#ifdef DCA_WITH_QMC_BIT
  //  pthread_mutex_lock(&mutex_numerical_error);
  // accumulator.get_error_distribution() += walker.get_error_distribution();
  //  pthread_mutex_unlock(&mutex_numerical_error);
  //#endif  // DCA_WITH_QMC_BIT

  if (id == 0 && concurrency.id() == concurrency.first()) {
    std::cout << "\n\t\t QMCI ends\n" << std::endl;
    walker.printSummary();
  }

  profiler_type::stop_threading(id);
}

template <class qmci_integrator_type>
void StdThreadQmciClusterSolver<qmci_integrator_type>::warm_up(walker_type& walker, int id) {
  if (id == 0) {
    if (concurrency.id() == concurrency.first())
      std::cout << "\n\t\t warm-up starts\n" << std::endl;
  }

  for (int i = 0; i < parameters.get_warm_up_sweeps(); i++) {
    walker.do_sweep();

    if (id == 0)
      walker.update_shell(i, parameters.get_warm_up_sweeps());
  }

  walker.is_thermalized() = true;

  if (id == 0) {
    if (concurrency.id() == concurrency.first())
      std::cout << "\n\t\t warm-up ends\n" << std::endl;
  }
}

template <class qmci_integrator_type>
void StdThreadQmciClusterSolver<qmci_integrator_type>::start_accumulator(int id) {
  profiler_type::start_threading(id);

  const int n_meas =
      parallel::util::getWorkload(parameters.get_measurements(), parameters.get_accumulators(),
                                  thread_task_handler_.IDToAccumIndex(id), concurrency);

  stdthread_accumulator_type accumulator_obj(parameters, data_, n_meas, id);

  accumulator_obj.initialize(DCA_iteration);

  for (int i = 0; i < n_meas; ++i) {
    {
      std::lock_guard<std::mutex> lock(mutex_queue);
      accumulators_queue.push(&accumulator_obj);
    }
    queue_insertion_.notify_one();

    {
      profiler_type profiler("stdthread-accumulator waiting", "stdthread-MC-accumulator", __LINE__,
                             id);
      accumulator_obj.wait_for_qmci_walker();
    }

    {
      profiler_type profiler("stdthread-accumulator accumulating", "stdthread-MC-accumulator",
                             __LINE__, id);
      accumulator_obj.measure();
    }
  }

  ++acc_finished;
  {
    std::lock_guard<std::mutex> lock(mutex_merge);
    accumulator_obj.sum_to(accumulator);
  }

  profiler_type::stop_threading(id);
}

template <class qmci_integrator_type>
void StdThreadQmciClusterSolver<qmci_integrator_type>::start_walker_and_accumulator(int id) {
  profiler_type::start_threading(id);

  // Create and warm a walker.
  walker_type walker(parameters, data_, rng_vector[id], id);
  walker.initialize();
  {
    profiler_type profiler("thermalization", "stdthread-MC", __LINE__, id);
    warm_up(walker, id);
  }

  accumulator_type accumulator_obj(parameters, data_, id);
  accumulator_obj.initialize(DCA_iteration);

  const int n_meas = parallel::util::getWorkload(parameters.get_measurements(),
                                                 parameters.get_accumulators(), id, concurrency);

  for (int i = 0; i < n_meas; ++i) {
    {
      profiler_type profiler("Walker updating", "stdthread-MC", __LINE__, id);
      walker.do_sweep();
    }
    {
      profiler_type profiler("Accumulator measuring", "stdthread-MC", __LINE__, id);
      accumulator_obj.update_from(walker);
      accumulator_obj.measure();
    }
    if (id == 0)
      walker.update_shell(i, n_meas);
  }

  ++acc_finished;
  std::lock_guard<std::mutex> lock(mutex_merge);
  accumulator_obj.sum_to(accumulator);

  profiler_type::stop_threading(id);
}

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_CLUSTER_SOLVER_HPP
