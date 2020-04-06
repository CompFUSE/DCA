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
#include <future>
#include <queue>
#include <stdexcept>
#include <vector>

#include "dca/io/buffer.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/linalg/util/handle_functions.hpp"
#include "dca/parallel/stdthread/thread_pool/thread_pool.hpp"
#include "dca/parallel/util/get_workload.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/qmci_autocorrelation_data.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_walker.hpp"
#include "dca/phys/dca_step/cluster_solver/thread_task_handler.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <class QmciSolver>
class StdThreadQmciClusterSolver : public QmciSolver {
public:
  using BaseClass = QmciSolver;
  using ThisType = StdThreadQmciClusterSolver<BaseClass>;

  using Data = typename BaseClass::DataType;
  using Parameters = typename BaseClass::ParametersType;
  using typename BaseClass::Concurrency;
  using typename BaseClass::Profiler;
  using typename BaseClass::Rng;

  using typename BaseClass::Accumulator;
  using Walker = stdthreadqmci::StdThreadQmciWalker<typename BaseClass::Walker>;
  using StdThreadAccumulatorType = stdthreadqmci::StdThreadQmciAccumulator<Accumulator>;

  StdThreadQmciClusterSolver(Parameters& parameters_ref, Data& data_ref,
                             io::HDF5Writer* file = nullptr);

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

private:
  void startWalker(int id);
  void startAccumulator(int id);
  void startWalkerAndAccumulator(int id);

  void initializeAndWarmUp(Walker& walker, int id, int walker_id);

  void readConfigurations();
  void writeConfigurations() const;

  void iterateOverLocalMeasurements(int walker_id, std::function<void(int, int, bool)>&& f);

  void printIntegrationMetadata() const;

private:
  using BaseClass::accumulator_;
  using BaseClass::concurrency_;
  using BaseClass::data_;
  using BaseClass::dca_iteration_;
  using BaseClass::parameters_;
  using BaseClass::total_time_;

  std::atomic<int> walk_finished_;
  std::atomic<uint> measurements_done_;

  const int nr_walkers_;
  const int nr_accumulators_;
  std::vector<std::size_t> walker_fingerprints_;
  std::vector<std::size_t> accum_fingerprints_;

  ThreadTaskHandler thread_task_handler_;

  std::vector<Rng> rng_vector_;

  std::queue<StdThreadAccumulatorType*> accumulators_queue_;

  std::mutex mutex_merge_;
  std::mutex mutex_queue_;
  std::condition_variable queue_insertion_;

  std::vector<dca::io::Buffer> config_dump_;
  stdthreadqmci::QmciAutocorrelationData<typename BaseClass::Walker> autocorrelation_data_;

  io::HDF5Writer* writer_ = nullptr;

  bool last_iteration_ = false;
  int measurements_ = 0;
};

template <class QmciSolver>
StdThreadQmciClusterSolver<QmciSolver>::StdThreadQmciClusterSolver(Parameters& parameters_ref,
                                                                   Data& data_ref,
                                                                   io::HDF5Writer* writer)
    : BaseClass(parameters_ref, data_ref),

      nr_walkers_(parameters_.get_walkers()),
      nr_accumulators_(parameters_.get_accumulators()),

      walker_fingerprints_(nr_walkers_, 0),
      accum_fingerprints_(nr_accumulators_, 0),

      thread_task_handler_(nr_walkers_, nr_accumulators_,
                           parameters_ref.shared_walk_and_accumulation_thread()),

      accumulators_queue_(),

      config_dump_(nr_walkers_),
      autocorrelation_data_(parameters_, 0),

      writer_(writer) {
  if (nr_walkers_ < 1 || nr_accumulators_ < 1) {
    throw std::logic_error(
        "Both the number of walkers and the number of accumulators must be at least 1.");
  }

  for (int i = 0; i < nr_walkers_; ++i) {
    rng_vector_.emplace_back(concurrency_.id(), concurrency_.number_of_processors(),
                             parameters_.get_seed());
  }

  readConfigurations();

  // Create a sufficient amount of cublas handles, cuda streams and threads.
  linalg::util::resizeHandleContainer(thread_task_handler_.size());
  parallel::ThreadPool::get_instance().enlarge(thread_task_handler_.size());
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::initialize(int dca_iteration) {
  Profiler profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);

  last_iteration_ = dca_iteration == parameters_.get_dca_iterations() - 1;

  measurements_ = parameters_.get_measurements()[dca_iteration];

  BaseClass::initialize(dca_iteration);

  walk_finished_ = 0;
  measurements_done_ = 0;
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::integrate() {
  Profiler profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);

  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "Threaded QMC integration has started: " << dca::util::print_time() << "\n"
              << std::endl;
  }

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

  auto print_metadata = [&]() {
    assert(walk_finished_ == parameters_.get_walkers());

    dca::profiling::WallTime end_time;

    dca::profiling::Duration duration(end_time, start_time);
    total_time_ = duration.sec + 1.e-6 * duration.usec;

    printIntegrationMetadata();
  };

  try {
    for (auto& future : futures)
      future.get();
  }
  catch (std::exception& err) {
    print_metadata();
    throw;
  }

  print_metadata();

  QmciSolver::accumulator_.finalize();
}

template <class QmciSolver>
template <typename dca_info_struct_t>
double StdThreadQmciClusterSolver<QmciSolver>::finalize(dca_info_struct_t& dca_info_struct) {
  Profiler profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);
  if (dca_iteration_ == parameters_.get_dca_iterations() - 1)
    BaseClass::computeErrorBars();

  double L2_Sigma_difference = QmciSolver::finalize(dca_info_struct);

  if (dca_iteration_ == parameters_.get_dca_iterations() - 1)
    writeConfigurations();

  // Write and reset autocorrelation.
  autocorrelation_data_.sumConcurrency(concurrency_);
  if (writer_)
    autocorrelation_data_.write(*writer_, dca_iteration_);
  autocorrelation_data_.reset();

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

  auto walker_log = last_iteration_ ? writer_ : nullptr;
  Walker walker(parameters_, data_, rng_vector_[walker_index], id, walker_log);

  std::unique_ptr<std::exception> exception_ptr;

  try {
    initializeAndWarmUp(walker, id, walker_index);

    iterateOverLocalMeasurements(
        walker_index, [&](const int meas_id, const int tot_meas, const bool print) {
          StdThreadAccumulatorType* acc_ptr = nullptr;

          {
            Profiler profiler("stdthread-MC-walker updating", "stdthread-MC-walker", __LINE__, id);
            walker.doSweep();
          }
          if (print)
            walker.updateShell(meas_id, tot_meas);

          {
            Profiler profiler("stdthread-MC-walker waiting", "stdthread-MC-walker", __LINE__, id);
            acc_ptr = nullptr;

            // Wait for available accumulators.
            {
              std::unique_lock<std::mutex> lock(mutex_queue_);
              queue_insertion_.wait(lock, [&]() { return !accumulators_queue_.empty(); });
              acc_ptr = accumulators_queue_.front();
              accumulators_queue_.pop();
            }
          }
          acc_ptr->updateFrom(walker);
        });
  }
  catch (std::bad_alloc& err) {
    --measurements_done_;
    std::cerr << "Walker " << id << " crashed."
              << "Walker fingerprint: " << walker.deviceFingerprint() * 1e-6 << " MB." << std::endl;
    // The exception is fatal if no other walker can take on the work of this one.
    if (parameters_.fix_meas_per_walker() || walk_finished_ == parameters_.get_walkers() - 1)
      exception_ptr = std::make_unique<std::bad_alloc>(err);
  }

  // If this is the last walker signal to all the accumulators to exit the loop.
  if (++walk_finished_ == parameters_.get_walkers()) {
    std::lock_guard<std::mutex> lock(mutex_queue_);
    while (!accumulators_queue_.empty()) {
      accumulators_queue_.front()->notifyDone();
      accumulators_queue_.pop();
    }
  }

  if (id == 0 && concurrency_.id() == concurrency_.first()) {
    std::cout << "\n\t\t QMCI ends\n" << std::endl;
    walker.printSummary();
  }

  if (parameters_.store_configuration() || (parameters_.get_directory_config_write() != "" &&
                                            dca_iteration_ == parameters_.get_dca_iterations() - 1))
    config_dump_[walker_index] = walker.dumpConfig();

  walker_fingerprints_[walker_index] = walker.deviceFingerprint();

  autocorrelation_data_ += walker;

  Profiler::stop_threading(id);

  if (exception_ptr)
    throw(*exception_ptr);
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::initializeAndWarmUp(Walker& walker, int id,
                                                                 int walker_id) {
  Profiler profiler("thermalization", "stdthread-MC-walker", __LINE__, id);

  // Read previous configuration.
  if (config_dump_[walker_id].size()) {
    walker.readConfig(config_dump_[walker_id]);
    config_dump_[walker_id].setg(0);  // Ready to read again if it is not overwritten.
  }

  walker.initialize(dca_iteration_);

  if (id == 0 && concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t warm-up starts\n" << std::endl;

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

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::iterateOverLocalMeasurements(
    const int walker_id, std::function<void(int, int, bool)>&& f) {
  const bool fix_thread_meas = parameters_.fix_meas_per_walker();
  const int total_meas = parallel::util::getWorkload(measurements_, concurrency_);

  const int n_local_meas =
      fix_thread_meas ? parallel::util::getWorkload(total_meas, parameters_.get_walkers(), walker_id)
                      : total_meas;
  const bool print = fix_thread_meas ? walker_id == 0 : true;

  if (fix_thread_meas) {
    // Perform a fixed amount of loops with a private counter.
    for (int meas_id = 0; meas_id < n_local_meas; ++meas_id)
      f(meas_id, n_local_meas, print);
  }
  else {
    // Perform the total number of loop with a shared atomic counter.
    for (int meas_id = measurements_done_++; meas_id < n_local_meas; meas_id = measurements_done_++)
      f(meas_id, n_local_meas, print);
  }
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::startAccumulator(int id) {
  Profiler::start_threading(id);

  StdThreadAccumulatorType accumulator_obj(parameters_, data_, id);

  accumulator_obj.initialize(dca_iteration_);

  std::unique_ptr<std::exception> exception_ptr;

  try {
    while (true) {
      {
        std::lock_guard<std::mutex> lock(mutex_queue_);
        if (walk_finished_ == parameters_.get_walkers())
          break;
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
  }
  catch (std::bad_alloc& err) {
    --measurements_done_;
    std::cerr << "Accumulator " << id << " crashed.\n"
              << "Accumulator fingerprint: " << accumulator_obj.deviceFingerprint() * 1e-6 << " MB."
              << std::endl;
    if (parameters_.fix_meas_per_walker() || walk_finished_ == parameters_.get_walkers())
      exception_ptr = std::make_unique<std::bad_alloc>(err);
  }

  {
    std::lock_guard<std::mutex> lock(mutex_merge_);
    accumulator_obj.sumTo(QmciSolver::accumulator_);
  }

  accum_fingerprints_[thread_task_handler_.IDToAccumIndex(id)] = accumulator_obj.deviceFingerprint();
  Profiler::stop_threading(id);

  if (exception_ptr)
    throw(*exception_ptr);
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::startWalkerAndAccumulator(int id) {
  Profiler::start_threading(id);

  // Create and warm a walker.
  auto walker_log = last_iteration_ ? writer_ : nullptr;
  Walker walker(parameters_, data_, rng_vector_[id], id, walker_log);
  initializeAndWarmUp(walker, id, id);

  Accumulator accumulator_obj(parameters_, data_, id);
  accumulator_obj.initialize(dca_iteration_);

  std::unique_ptr<std::exception> current_exception;

  try {
    iterateOverLocalMeasurements(id, [&](const int meas_id, const int n_meas, const bool print) {
      {
        Profiler profiler("Walker updating", "stdthread-MC", __LINE__, id);
        walker.doSweep();
      }
      {
        Profiler profiler("Accumulator measuring", "stdthread-MC", __LINE__, id);
        accumulator_obj.updateFrom(walker);
        accumulator_obj.measure();
      }
      if (print)
        walker.updateShell(meas_id, n_meas);
    });
  }
  catch (std::bad_alloc& err) {
    --measurements_done_;
    std::cerr << "Walker and accumulator " << id << " crashed.\n"
              << "Walker fingerprint: " << walker.deviceFingerprint() * 1e-6 << " MB.\n"
              << "Accumulator fingerprint: " << accumulator_obj.deviceFingerprint() * 1e-6 << " MB."
              << std::endl;

    // The exception is fatal if no other walker can take on the work of this one.
    if (parameters_.fix_meas_per_walker() || walk_finished_ == parameters_.get_walkers() - 1)
      current_exception = std::make_unique<std::bad_alloc>(err);
  }

  ++walk_finished_;
  {
    std::lock_guard<std::mutex> lock(mutex_merge_);
    accumulator_obj.sumTo(QmciSolver::accumulator_);
  }

  if (parameters_.store_configuration() || (parameters_.get_directory_config_write() != "" &&
                                            dca_iteration_ == parameters_.get_dca_iterations() - 1))
    config_dump_[id] = walker.dumpConfig();

  walker_fingerprints_[id] = walker.deviceFingerprint();
  accum_fingerprints_[id] = accumulator_obj.deviceFingerprint();

  autocorrelation_data_ += walker;

  Profiler::stop_threading(id);

  if (current_exception)
    throw(*current_exception);
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::writeConfigurations() const {
  if (parameters_.get_directory_config_write() == "")
    return;

  try {
    const std::string out_name = parameters_.get_directory_config_write() + "/process_" +
                                 std::to_string(concurrency_.id()) + ".hdf5";
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
  if (parameters_.get_directory_config_read() == "")
    return;

  Profiler profiler(__FUNCTION__, "stdthread-MC", __LINE__);

  try {
    const std::string inp_name = parameters_.get_directory_config_read() + "/process_" +
                                 std::to_string(concurrency_.id()) + ".hdf5";
    io::HDF5Reader reader(false);
    reader.open_file(inp_name);
    for (int id = 0; id < config_dump_.size(); ++id)
      reader.execute("configuration_" + std::to_string(id), config_dump_[id]);
  }
  catch (std::exception& err) {
    std::cerr << err.what() << "\nCould not read the configuration.\n";
    for (auto& config : config_dump_)
      config.clear();
  }
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::printIntegrationMetadata() const {
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "Threaded on-node integration has ended: " << dca::util::print_time()
              << "\n\nTotal number of measurements: " << measurements_ << "\nQMC-time\t"
              << total_time_ << "\n";
    if (QmciSolver::device == linalg::GPU) {
      std::cout << "\nWalker fingerprints [MB]: \n";
      for (const auto& x : walker_fingerprints_)
        std::cout << x * 1e-6 << "\n";
      std::cout << "Accumulator fingerprints [MB]: \n";
      for (const auto& x : accum_fingerprints_)
        std::cout << x * 1e-6 << "\n";
      std::cout << "Static Accumulator fingerprint [MB]:\n"
                << Accumulator::staticDeviceFingerprint() * 1e-6 << "\n\n";
    }
  }
}

}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_CLUSTER_SOLVER_HPP
