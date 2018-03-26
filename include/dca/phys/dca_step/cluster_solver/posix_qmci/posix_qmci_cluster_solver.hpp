// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// A posix MC integrator that implements a threaded MC integration independent of the MC method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_POSIX_QMCI_CLUSTER_SOLVER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_POSIX_QMCI_CLUSTER_SOLVER_HPP

#include <pthread.h>

#include <iostream>
#include <queue>
#include <stdexcept>
#include <vector>

#include "dca/parallel/pthreading/lock.hpp"
#include "dca/parallel/pthreading/conditional_variable.hpp"
#include "dca/parallel/pthreading/pthreading.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_solver/posix_qmci/posix_qmci_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/posix_qmci/thread_task_handler.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <class QmcSolver>
class PosixQmciClusterSolver : public QmcSolver {
  using BaseClass = QmcSolver;
  using this_type = PosixQmciClusterSolver<BaseClass>;
  using typename BaseClass::Data;
  using typename BaseClass::Concurrency;
  using typename BaseClass::Parameters;
  using typename BaseClass::Profiler;
  using typename BaseClass::Rng;

public:
  PosixQmciClusterSolver(Parameters& parameters_ref, Data& data_ref);

  void initialize(int dca_iteration = 0);

  void integrate();

  using BaseClass::finalize;
  using BaseClass::write;

private:
  using typename BaseClass::Walker;
  using typename BaseClass::Accumulator;
  using BaseClass::printStatus;
  using PosixAccumulator = posixqmci::PosixQmciAccumulator<Accumulator>;

private:
  static void* startWalkerStatic(void* arg);
  static void* startAccumulatorStatic(void* arg);
  static void* startWalkerAndAccumulatorStatic(void* arg);

  void startWalker(int id);
  void startAccumulator(int id);
  void startWalkerAndAccumulator(int id);

  void sumThreadResults(Accumulator& accumulator);
  void warmUp(Walker& walker, int id);

private:
  using BaseClass::parameters_;
  using BaseClass::data_;
  using BaseClass::concurrency_;
  using BaseClass::total_time_;
  using BaseClass::dca_iteration_;

  using BaseClass::accumulator_;

  const int nr_walkers_;
  const int nr_accumulators_;
  bool shared_walk_acc_thread_;

  posixqmci::ThreadTaskHandler thread_task_handler_;
  std::vector<Rng> rng_vector_;

  std::queue<PosixAccumulator*> accumulators_queue_;

  ushort acc_finished_;

  parallel::Lock<parallel::Pthreading> merge_lock_;
  parallel::Lock<parallel::Pthreading> queue_lock_;
};

template <class QmcSolver>
PosixQmciClusterSolver<QmcSolver>::PosixQmciClusterSolver(Parameters& parameters_ref, Data& data_)
    : BaseClass(parameters_ref, data_),

      nr_walkers_(parameters_.get_walkers()),
      nr_accumulators_(parameters_.get_accumulators()),
      shared_walk_acc_thread_(parameters_.shared_walk_and_accumulation_thread()),
      thread_task_handler_(nr_walkers_, nr_accumulators_, shared_walk_acc_thread_),
      acc_finished_(0) {
  if (nr_walkers_ < 1 || nr_accumulators_ < 1) {
    throw std::logic_error(
        "Both the number of walkers and the number of accumulators must be at least 1.");
  }
  if ((nr_walkers_ != nr_accumulators_) and shared_walk_acc_thread_) {
    if (concurrency_.id() == concurrency_.first())
      std::cerr << "To use shared thread functionality the number of walkers and accumulators must "
                   "be the same. Disabling shared thread"
                << std::endl;
    shared_walk_acc_thread_ = false;
  }

  for (int i = 0; i < nr_walkers_; ++i) {
    rng_vector_.emplace_back(concurrency_.id(), concurrency_.number_of_processors(),
                             parameters_.get_seed());
  }
}

template <class QmcSolver>
void PosixQmciClusterSolver<QmcSolver>::initialize(int dca_iteration) {
  Profiler profiler(__FUNCTION__, "posix-MC-Integration", __LINE__);

  BaseClass::initialize(dca_iteration, false);

  accumulator_.release();
  acc_finished_ = 0;
}

template <class QmcSolver>
void PosixQmciClusterSolver<QmcSolver>::integrate() {
  Profiler profiler(__FUNCTION__, "posix-MC-Integration", __LINE__);

  if (concurrency_.id() == concurrency_.first())
    std::cout << "Threaded QMC integration has started: " << dca::util::print_time() << "\n"
              << std::endl;

  const int num_threads = shared_walk_acc_thread_ ? nr_walkers_ : nr_accumulators_ + nr_walkers_;
  std::vector<pthread_t> threads(num_threads);
  std::vector<std::pair<this_type*, int>> data(num_threads);

  if (concurrency_.id() == concurrency_.first())
    thread_task_handler_.print();

  dca::profiling::WallTime start_time;

  for (int i = 0; i < num_threads; ++i) {
    data[i] = std::pair<this_type*, int>(this, i);

    if (thread_task_handler_.getTask(i) == "walker")
      pthread_create(&threads[i], nullptr, startWalkerStatic, &data[i]);
    else if (thread_task_handler_.getTask(i) == "accumulator")
      pthread_create(&threads[i], nullptr, startAccumulatorStatic, &data[i]);
    else if (thread_task_handler_.getTask(i) == "walker and accumulator")
      pthread_create(&threads[i], nullptr, startWalkerAndAccumulatorStatic, &data[i]);
    else
      throw std::logic_error("Thread task is undefined.");
  }

  void* rc;
  for (int i = 0; i < num_threads; ++i)
    pthread_join(threads[i], &rc);

  dca::profiling::WallTime end_time;

  dca::profiling::Duration duration(end_time, start_time);
  total_time_ = duration.sec + 1.e-6 * duration.usec;

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t threaded QMC integration ends\n" << std::endl;

  if (concurrency_.id() == concurrency_.first())
    std::cout << "Threaded on-node integration has ended: " << dca::util::print_time()
              << "\n\nTotal number of measurements: "
              << concurrency_.number_of_processors() *
                     parameters_.get_measurements_per_process_and_accumulator() * nr_accumulators_
              << std::endl;

  accumulator_->finalize();
}

template <class QmcSolver>
void* PosixQmciClusterSolver<QmcSolver>::startWalkerStatic(void* arg) {
  std::pair<this_type*, int>* data = reinterpret_cast<std::pair<this_type*, int>*>(arg);

  Profiler::start_pthreading(data->second);

  data->first->startWalker(data->second);

  Profiler::stop_pthreading(data->second);

  return NULL;
}

template <class QmcSolver>
void* PosixQmciClusterSolver<QmcSolver>::startAccumulatorStatic(void* arg) {
  std::pair<this_type*, int>* data = reinterpret_cast<std::pair<this_type*, int>*>(arg);

  Profiler::start_pthreading(data->second);

  data->first->startAccumulator(data->second);

  Profiler::stop_pthreading(data->second);

  return NULL;
}

template <class QmcSolver>
void* PosixQmciClusterSolver<QmcSolver>::startWalkerAndAccumulatorStatic(void* arg) {
  std::pair<this_type*, int>* data = reinterpret_cast<std::pair<this_type*, int>*>(arg);
  Profiler::start_pthreading(data->second);
  data->first->startWalkerAndAccumulator(data->second);
  Profiler::stop_pthreading(data->second);
  return NULL;
}

template <class QmcSolver>
void PosixQmciClusterSolver<QmcSolver>::startWalkerAndAccumulator(int id) {
  if (id == 0 and concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t QMCI starts\n" << std::endl;

  const int rng_index = thread_task_handler_.walkerIDToRngIndex(id);
  Walker walker(BaseClass::instantiateWalker(rng_vector_[rng_index], id));
  {
    Profiler profiler("Thermalization", "posix-MC-walker", __LINE__, id);
    warmUp(walker, id);
  }
  Accumulator accumulator(BaseClass::instantiateAccumulator(id));

  const int n_meas = parameters_.get_measurements_per_process_and_accumulator();

  for (int i = 0; i < n_meas; ++i) {
    {
      Profiler profiler("Walker updating", "posix-MC-walker", __LINE__, id);
      walker.doSweep();
    }
    {
      Profiler profiler("Accumulator measuring", "posix-MC-accumulator", __LINE__, id);
      accumulator.accumulate(walker);
    }
    if (id == 0)
      printStatus(i, n_meas, walker);
  }

  sumThreadResults(accumulator);
}

template <class QmcSolver>
void PosixQmciClusterSolver<QmcSolver>::sumThreadResults(Accumulator& accumulator) {
  merge_lock_.lock();
  if (not acc_finished_)
    accumulator_.reset(new Accumulator(std::move(accumulator)));
  else
    accumulator.sumTo(*accumulator_);
  ++acc_finished_;
  merge_lock_.unlock();
}

template <class QmcSolver>
void PosixQmciClusterSolver<QmcSolver>::warmUp(Walker& walker, int id) {
  if (id == 0 and concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t warm-up starts\n" << std::endl;

  const int meas_to_do = parameters_.get_warm_up_sweeps();
  for (int i = 0; i < meas_to_do; i++) {
    walker.doSweep();

    if (id == 0)
      printStatus(i, meas_to_do, walker);
  }

  walker.markThermalized();

  if (id == 0 and concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t warm-up ends\n" << std::endl;
}

template <class QmcSolver>
void PosixQmciClusterSolver<QmcSolver>::startWalker(int id) {
  if (id == 0 and concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t QMCI starts\n\n";

  const int rng_index = thread_task_handler_.walkerIDToRngIndex(id);
  Walker walker(BaseClass::instantiateWalker(rng_vector_[rng_index], id));
  {
    Profiler profiler("Thermalization", "posix-MC-walker", __LINE__, id);
    warmUp(walker, id);
  }

  PosixAccumulator* acc_ptr(nullptr);

  while (acc_finished_ < nr_accumulators_) {
    {
      Profiler profiler("posix-MC-walker updating", "posix-MC-walker", __LINE__, id);
      walker.doSweep();
    }

    {
      Profiler profiler("posix-MC-walker waiting", "posix-MC-walker", __LINE__, id);

      while (acc_finished_ < nr_accumulators_) {
        acc_ptr = nullptr;

        // checking for available accumulators
        queue_lock_.lock();
        if (!accumulators_queue_.empty()) {
          acc_ptr = accumulators_queue_.front();
          accumulators_queue_.pop();
        }
        queue_lock_.unlock();

        if (acc_ptr != nullptr) {
          acc_ptr->updateFrom(walker);
          acc_ptr = nullptr;
          break;
        }
      }
    }
  }
}

template <class QmcSolver>
void PosixQmciClusterSolver<QmcSolver>::startAccumulator(int id) {
  PosixAccumulator accumulator_obj(BaseClass::instantiateAccumulator(id));
  const int meas_to_do = parameters_.get_measurements_per_process_and_accumulator();

  for (int i = 0; i < meas_to_do; ++i) {
    queue_lock_.lock();
    accumulators_queue_.push(&accumulator_obj);
    queue_lock_.unlock();
    {
      Profiler profiler("posix-accumulator waiting", "posix-MC-accumulator", __LINE__, id);
      accumulator_obj.waitForQmciWalker();
    }
    {
      Profiler profiler("posix-accumulator accumulating", "posix-MC-accumulator", __LINE__, id);
      accumulator_obj.measure();
      // Update shell.
      if (id == 0)
        printStatus(accumulator_obj.get_number_of_measurements(), meas_to_do,
                    static_cast<Accumulator&>(accumulator_obj));
    }
  }

  sumThreadResults(static_cast<Accumulator&>(accumulator_obj));
}

}  // solver
}  // phys
}  // dca

//#endif  // DCA_HAVE_PTHREADS
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_POSIX_CTINT_CLUSTER_SOLVER_HPP
