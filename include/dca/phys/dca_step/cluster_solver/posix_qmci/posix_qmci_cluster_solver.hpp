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
//
// A posix MC integrator that implements a threaded MC integration independent of the MC method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_POSIX_QMCI_CLUSTER_SOLVER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_POSIX_QMCI_CLUSTER_SOLVER_HPP

#include <pthread.h>

#include <iostream>
#include <queue>
#include <stdexcept>
#include <vector>

#include "dca/phys/dca_step/cluster_solver/posix_qmci/posix_qmci_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/posix_qmci/thread_task_handler.hpp"
#include "dca/profiling/events/time.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <class qmci_integrator_type>
class PosixQmciClusterSolver : public qmci_integrator_type {
  typedef typename qmci_integrator_type::this_MOMS_type MOMS_type;
  typedef typename qmci_integrator_type::this_parameters_type parameters_type;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  using random_number_generator = typename parameters_type::random_number_generator;

  typedef typename qmci_integrator_type::walker_type walker_type;
  typedef typename qmci_integrator_type::accumulator_type accumulator_type;

  typedef PosixQmciClusterSolver<qmci_integrator_type> this_type;
  typedef posixqmci::posix_qmci_accumulator<accumulator_type> posix_accumulator_type;

public:
  PosixQmciClusterSolver(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  template <typename Writer>
  void write(Writer& reader);

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

private:
  static void* start_walker_static(void* arg);
  static void* start_accumulator_static(void* arg);

  void start_walker(int id);
  void start_accumulator(int id);

  void warm_up(walker_type& walker, int id);

  // TODO: Are the following using statements redundant and can therefore be removed?
  using qmci_integrator_type::compute_error_bars;
  using qmci_integrator_type::symmetrize_measurements;

private:
  using qmci_integrator_type::parameters;
  using qmci_integrator_type::MOMS;
  using qmci_integrator_type::concurrency;

  using qmci_integrator_type::total_time;

  using qmci_integrator_type::DCA_iteration;

  using qmci_integrator_type::accumulator;

  int acc_finished;

  const int nr_walkers;
  const int nr_accumulators;

  posixqmci::ThreadTaskHandler thread_task_handler_;

  std::vector<random_number_generator> rng_vector;

  std::queue<posix_accumulator_type*> accumulators_queue;

  pthread_mutex_t mutex_print;
  pthread_mutex_t mutex_merge;
  pthread_mutex_t mutex_queue;
  pthread_mutex_t mutex_acc_finished;

  pthread_mutex_t mutex_numerical_error;
};

template <class qmci_integrator_type>
PosixQmciClusterSolver<qmci_integrator_type>::PosixQmciClusterSolver(parameters_type& parameters_ref,
                                                                     MOMS_type& MOMS_ref)
    : qmci_integrator_type(parameters_ref, MOMS_ref),

      nr_walkers(parameters.get_walkers()),
      nr_accumulators(parameters.get_accumulators()),

      thread_task_handler_(nr_walkers, nr_accumulators),

      accumulators_queue() {
  if (nr_walkers < 1 || nr_accumulators < 1) {
    throw std::logic_error(
        "Both the number of walkers and the number of accumulators must be at least 1.");
  }

  for (int i = 0; i < nr_walkers; ++i) {
    rng_vector.emplace_back(concurrency.id(), concurrency.number_of_processors(),
                            parameters.get_seed());
  }
}

template <class qmci_integrator_type>
template <typename Writer>
void PosixQmciClusterSolver<qmci_integrator_type>::write(Writer& writer) {
  qmci_integrator_type::write(writer);
  // accumulator.write(writer);
}

template <class qmci_integrator_type>
void PosixQmciClusterSolver<qmci_integrator_type>::initialize(int dca_iteration) {
  profiler_type profiler(__FUNCTION__, "posix-MC-Integration", __LINE__);

  qmci_integrator_type::initialize(dca_iteration);

  acc_finished = 0;

  pthread_mutex_init(&mutex_print, NULL);
  pthread_mutex_init(&mutex_merge, NULL);
  pthread_mutex_init(&mutex_queue, NULL);

  pthread_mutex_init(&mutex_acc_finished, NULL);
  pthread_mutex_init(&mutex_numerical_error, NULL);
}

template <class qmci_integrator_type>
void PosixQmciClusterSolver<qmci_integrator_type>::integrate() {
  profiler_type profiler(__FUNCTION__, "posix-MC-Integration", __LINE__);

  if (concurrency.id() == concurrency.first()) {
    std::cout << "Threaded QMC integration has started: " << dca::util::print_time() << "\n"
              << std::endl;
  }

  std::vector<pthread_t> threads(nr_accumulators + nr_walkers);
  std::vector<std::pair<this_type*, int>> data(nr_accumulators + nr_walkers);

  {
    if (concurrency.id() == concurrency.first())
      thread_task_handler_.print();

    dca::profiling::WallTime start_time;

    for (int i = 0; i < nr_walkers + nr_accumulators; ++i) {
      data[i] = std::pair<this_type*, int>(this, i);

      if (thread_task_handler_.getTask(i) == "walker")
        pthread_create(&threads[i], NULL, start_walker_static, &data[i]);
      else if (thread_task_handler_.getTask(i) == "accumulator")
        pthread_create(&threads[i], NULL, start_accumulator_static, &data[i]);
      else
        throw std::logic_error("Thread is neither a walker nor an accumulator.");
    }

    void* rc;
    for (int i = 0; i < nr_walkers + nr_accumulators; ++i) {
      pthread_join(threads[i], &rc);
    }

    dca::profiling::WallTime end_time;

    dca::profiling::Duration duration(end_time, start_time);
    total_time = duration.sec + 1.e-6 * duration.usec;
  }

  if (concurrency.id() == concurrency.first()) {
    std::cout << "Threaded on-node integration has ended: " << dca::util::print_time()
              << "\n\nTotal number of measurements: "
              << concurrency.number_of_processors() *
                     parameters.get_measurements_per_process_and_accumulator() * nr_accumulators
              << std::endl;
  }
}

template <class qmci_integrator_type>
template <typename dca_info_struct_t>
double PosixQmciClusterSolver<qmci_integrator_type>::finalize(dca_info_struct_t& dca_info_struct) {
  profiler_type profiler(__FUNCTION__, "posix-MC-Integration", __LINE__);
  if (DCA_iteration == parameters.get_dca_iterations() - 1)
    compute_error_bars();

  double L2_Sigma_difference = qmci_integrator_type::finalize(dca_info_struct);

  pthread_mutex_destroy(&mutex_print);
  pthread_mutex_destroy(&mutex_merge);
  pthread_mutex_destroy(&mutex_queue);

  pthread_mutex_destroy(&mutex_acc_finished);
  pthread_mutex_destroy(&mutex_numerical_error);

  return L2_Sigma_difference;
}

template <class qmci_integrator_type>
void* PosixQmciClusterSolver<qmci_integrator_type>::start_walker_static(void* arg) {
  std::pair<this_type*, int>* data = reinterpret_cast<std::pair<this_type*, int>*>(arg);

  profiler_type::start_pthreading(data->second);

  data->first->start_walker(data->second);

  profiler_type::stop_pthreading(data->second);

  return NULL;
}

template <class qmci_integrator_type>
void* PosixQmciClusterSolver<qmci_integrator_type>::start_accumulator_static(void* arg) {
  std::pair<this_type*, int>* data = reinterpret_cast<std::pair<this_type*, int>*>(arg);

  profiler_type::start_pthreading(data->second);

  data->first->start_accumulator(data->second);

  profiler_type::stop_pthreading(data->second);

  return NULL;
}

template <class qmci_integrator_type>
void PosixQmciClusterSolver<qmci_integrator_type>::start_walker(int id) {
  if (id == 0) {
    if (concurrency.id() == concurrency.first())
      std::cout << "\n\t\t QMCI starts\n" << std::endl;
  }

  const int rng_index = thread_task_handler_.walkerIDToRngIndex(id);
  walker_type walker(parameters, MOMS, rng_vector[rng_index], id);

  walker.initialize();

  {
    profiler_type profiler("thermalization", "posix-MC-walker", __LINE__, id);
    warm_up(walker, id);
  }

  posix_accumulator_type* acc_ptr(NULL);

  while (acc_finished < nr_accumulators) {
    {
      profiler_type profiler("posix-MC-walker updating", "posix-MC-walker", __LINE__, id);
      walker.do_sweep();
    }

    {
      profiler_type profiler("posix-MC-walker waiting", "posix-MC-walker", __LINE__, id);

      while (acc_finished < nr_accumulators) {
        acc_ptr = NULL;

        {  // checking for available accumulators
          pthread_mutex_lock(&mutex_queue);

          if (!accumulators_queue.empty()) {
            acc_ptr = accumulators_queue.front();
            accumulators_queue.pop();
          }

          pthread_mutex_unlock(&mutex_queue);
        }

        if (acc_ptr != NULL) {
          acc_ptr->update_from(walker);
          acc_ptr = NULL;
          break;
        }

        for (int i = 0; i < parameters.get_additional_steps(); ++i) {
          // std::cout << "\twalker " << id << " is doing some additional steps !!\n";
          profiler_type profiler("additional steps", "posix-MC-walker", __LINE__, id);
          walker.do_step();
        }
      }
    }
  }

#ifdef DCA_WITH_QMC_BIT
  pthread_mutex_lock(&mutex_numerical_error);
  // accumulator.get_error_distribution() += walker.get_error_distribution();
  pthread_mutex_unlock(&mutex_numerical_error);
#endif  // DCA_WITH_QMC_BIT

  if (id == 0) {
    if (concurrency.id() == concurrency.first())
      std::cout << "\n\t\t QMCI ends\n" << std::endl;
  }
}

template <class qmci_integrator_type>
void PosixQmciClusterSolver<qmci_integrator_type>::warm_up(walker_type& walker, int id) {
  if (id == 0) {
    if (concurrency.id() == concurrency.first())
      std::cout << "\n\t\t warm-up starts\n" << std::endl;
  }

  for (int i = 0; i < parameters.get_warm_up_sweeps(); i++) {
    walker.do_sweep();

    if (id == 0)
      this->update_shell(i, parameters.get_warm_up_sweeps(), walker.get_configuration().size());
  }

  walker.is_thermalized() = true;

  if (id == 0) {
    if (concurrency.id() == concurrency.first())
      std::cout << "\n\t\t warm-up ends\n" << std::endl;
  }
}

template <class qmci_integrator_type>
void PosixQmciClusterSolver<qmci_integrator_type>::start_accumulator(int id) {
  posix_accumulator_type accumulator_obj(parameters, MOMS, id);

  accumulator_obj.initialize(DCA_iteration);

  for (int i = 0; i < parameters.get_measurements_per_process_and_accumulator(); ++i) {
    pthread_mutex_lock(&mutex_queue);
    accumulators_queue.push(&accumulator_obj);
    pthread_mutex_unlock(&mutex_queue);

    {
      profiler_type profiler("posix-accumulator waiting", "posix-MC-accumulator", __LINE__, id);
      accumulator_obj.wait_for_qmci_walker();
    }

    {
      profiler_type profiler("posix-accumulator accumulating", "posix-MC-accumulator", __LINE__, id);
      if (id == 1)
        this->update_shell(i, parameters.get_measurements_per_process_and_accumulator(),
                           accumulator_obj.get_configuration().size());

      accumulator_obj.measure(mutex_queue, accumulators_queue);
    }
  }

  {
    pthread_mutex_lock(&mutex_merge);

    acc_finished++;
    accumulator_obj.sum_to(accumulator);

    pthread_mutex_unlock(&mutex_merge);
  }
}

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_POSIX_QMCI_POSIX_QMCI_CLUSTER_SOLVER_HPP
