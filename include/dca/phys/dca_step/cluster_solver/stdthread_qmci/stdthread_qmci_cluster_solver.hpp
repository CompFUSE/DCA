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

#include "dca/io/buffer.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/thread_task_handler.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/util/get_stdout_from_command.hpp"
#include "dca/util/print_time.hpp"
#include "dca/parallel/util/get_workload.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <class qmci_integrator_type>
class StdThreadQmciClusterSolver : public qmci_integrator_type {
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
  static void* start_walker_static(pair_type data);
  static void* start_accumulator_static(pair_type data);

  void start_walker(int id);
  void start_accumulator(int id);

  void warm_up(walker_type& walker, int id);

  void readConfigurations();
  void writeConfigurations() const;
  int findAvailableFiles() const;

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
  std::atomic<int> acc_id_;
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

  std::vector<dca::io::Buffer> config_dump_;
};

template <class qmci_integrator_type>
StdThreadQmciClusterSolver<qmci_integrator_type>::StdThreadQmciClusterSolver(
    parameters_type& parameters_ref, Data& data_ref)
    : qmci_integrator_type(parameters_ref, data_ref),

      nr_walkers(parameters.get_walkers()),
      nr_accumulators(parameters.get_accumulators()),

      thread_task_handler_(nr_walkers, nr_accumulators),

      accumulators_queue(),
      config_dump_(nr_walkers) {
  if (nr_walkers < 1 || nr_accumulators < 1) {
    throw std::logic_error(
        "Both the number of walkers and the number of accumulators must be at least 1.");
  }

  for (int i = 0; i < nr_walkers; ++i) {
    rng_vector.emplace_back(concurrency.id(), concurrency.number_of_processors(),
                            parameters.get_seed());
  }

  readConfigurations();
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
  acc_id_ = 0;
}

template <class qmci_integrator_type>
void StdThreadQmciClusterSolver<qmci_integrator_type>::integrate() {
  profiler_type profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);

  if (concurrency.id() == concurrency.first()) {
    std::cout << "Threaded QMC integration has started: " << dca::util::print_time() << "\n"
              << std::endl;
  }

  measurements_remaining_ = parallel::util::getWorkload(parameters.get_measurements(), concurrency);
  ;

  std::vector<std::thread> threads;
  std::vector<pair_type> data;

  {
    if (concurrency.id() == concurrency.first())
      thread_task_handler_.print();

    dca::profiling::WallTime start_time;

    for (int i = 0; i < nr_walkers + nr_accumulators; ++i) {
      data.push_back(std::make_pair(this, i));

      if (thread_task_handler_.getTask(i) == "walker")
        threads.push_back(std::thread(start_walker_static, data.back()));
      else if (thread_task_handler_.getTask(i) == "accumulator")
        threads.push_back(std::thread(start_accumulator_static, data.back()));
      else
        throw std::logic_error("Thread is neither a walker nor an accumulator.");
    }

    for (int i = 0; i < nr_walkers + nr_accumulators; ++i) {
      threads[i].join();
    }

    dca::profiling::WallTime end_time;

    dca::profiling::Duration duration(end_time, start_time);
    total_time = duration.sec + 1.e-6 * duration.usec;
  }

  if (concurrency.id() == concurrency.first()) {
    std::cout << "Threaded on-node integration has ended: " << dca::util::print_time()
              << "\n\nTotal number of measurements: " << parameters.get_measurements() << std::endl;
  }
}

template <class qmci_integrator_type>
template <typename dca_info_struct_t>
double StdThreadQmciClusterSolver<qmci_integrator_type>::finalize(dca_info_struct_t& dca_info_struct) {
  profiler_type profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);
  if (DCA_iteration == parameters.get_dca_iterations() - 1)
    compute_error_bars();

  if (DCA_iteration == parameters.get_dca_iterations() - 1)
    writeConfigurations();

  double L2_Sigma_difference = qmci_integrator_type::finalize(dca_info_struct);
  return L2_Sigma_difference;
}

template <class qmci_integrator_type>
void* StdThreadQmciClusterSolver<qmci_integrator_type>::start_walker_static(pair_type data) {
  profiler_type::start_threading(data.second);

  data.first->start_walker(data.second);

  profiler_type::stop_threading(data.second);

  return NULL;
}

template <class qmci_integrator_type>
void* StdThreadQmciClusterSolver<qmci_integrator_type>::start_accumulator_static(pair_type data) {
  profiler_type::start_threading(data.second);

  data.first->start_accumulator(data.second);

  profiler_type::stop_threading(data.second);

  return NULL;
}

template <class qmci_integrator_type>
void StdThreadQmciClusterSolver<qmci_integrator_type>::start_walker(int id) {
  if (id == 0) {
    if (concurrency.id() == concurrency.first())
      std::cout << "\n\t\t QMCI starts\n" << std::endl;
  }

  const int rng_index = thread_task_handler_.walkerIDToRngIndex(id);
  walker_type walker(parameters, data_, rng_vector[rng_index], id);

  if (config_dump_[rng_index].size())
    walker.readConfig(config_dump_[rng_index]);

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
      acc_ptr = NULL;

      while (acc_ptr == NULL) {  // checking for available accumulators
        std::unique_lock<std::mutex> lock(mutex_queue);
        if (!accumulators_queue.empty()) {
          acc_ptr = accumulators_queue.front();
          accumulators_queue.pop();
        }
      }

      acc_ptr->update_from(walker);
    }
  }

#ifdef DCA_WITH_QMC_BIT
  mutex_numerical_error.lock();
  // accumulator.get_error_distribution() += walker.get_error_distribution();
  mutex_numerical_error.unlock();
#endif  // DCA_WITH_QMC_BIT

  config_dump_[rng_index] = walker.dumpConfig();

  if (id == 0 && concurrency.id() == concurrency.first()) {
    std::cout << "\n\t\t QMCI ends\n" << std::endl;
    walker.printSummary();
  }
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
  stdthread_accumulator_type accumulator_obj(parameters, data_, id);

  accumulator_obj.initialize(DCA_iteration);

  const int accumulator_id = acc_id_++;

  const int n_local_meas = parallel::util::getWorkload(
      parameters.get_measurements(), parameters.get_accumulators(), accumulator_id);

  for (int i = 0; i < n_local_meas; ++i) {
    {
      std::lock_guard<std::mutex> lock(mutex_queue);
      accumulators_queue.push(&accumulator_obj);
    }

    {
      profiler_type profiler("stdthread-accumulator waiting", "stdthread-MC-accumulator", __LINE__,
                             id);
      accumulator_obj.wait_for_qmci_walker();
    }

    {
      profiler_type profiler("stdthread-accumulator accumulating", "stdthread-MC-accumulator",
                             __LINE__, id);
      accumulator_obj.measure(mutex_queue, accumulators_queue);
    }
  }

  ++acc_finished;
  {
    std::lock_guard<std::mutex> lock(mutex_merge);
    accumulator_obj.sum_to(accumulator);
  }
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::writeConfigurations() const {
  if (parameters.get_directory_config_write() == "")
    return;

  try {
    const std::string out_name = parameters.get_directory_config_write() + "/process_" +
                                 std::to_string(concurrency.id()) + ".hdf5";
    io::HDF5Writer writer(false);
    writer.open_file(out_name);
    for (int id = 0; id < config_dump_.size(); ++id)
      writer.execute("configuration_" + std::to_string(id), config_dump_[id]);
  }
  catch (std::exception& err) {
    std::cerr << err.what() << "\nCould not write the configuration.\n";
  }
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::readConfigurations() {
  if (parameters.get_directory_config_read() == "")
    return;

  try {
    const int n_available = findAvailableFiles();
    const int id_to_read = concurrency.id() % n_available;

    const std::string inp_name = parameters.get_directory_config_read() + "/process_" +
                                 std::to_string(id_to_read) + ".hdf5";
    io::HDF5Reader reader(false);
    reader.open_file(inp_name);
    for (int id = 0; id < config_dump_.size(); ++id)
      reader.execute("configuration_" + std::to_string(id), config_dump_[id]);

    if (concurrency.id() == 0) {
      std::cout << "Read configuration from " << parameters.get_directory_config_read() << ".\n";
    }
  }
  catch (std::exception& err) {
    std::cerr << err.what() << "\nCould not read the configuration.\n";
    for (auto& config : config_dump_)
      config.clear();
  }
}

template <class QmciSolver>
int StdThreadQmciClusterSolver<QmciSolver>::findAvailableFiles() const {
  int result = 0;
  if (concurrency.id() == 0) {
    try {
      // Count the number of configuration files.
      const std::string cmd =
          "ls -1 " + parameters.get_directory_config_read() + "/process_*.hdf5 | wc -l";
      result = std::atoi(dca::util::getStdoutFromCommand(cmd).c_str());
    }
    catch (...) {
    }
  }

  concurrency.broadcast(result, 0);
  return result;
}

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_CLUSTER_SOLVER_HPP
