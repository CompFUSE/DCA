// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: John Biddiscombe (john.biddiscombe@cscs.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// An MC integrator that implements a threaded MC integration independent of the MC
// method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_CLUSTER_SOLVER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_CLUSTER_SOLVER_HPP

#include <atomic>
#include <iostream>
#include <future>
#include <queue>
#include <stdexcept>
#include <vector>

#include "dca/config/threading.hpp"
#include "dca/io/buffer.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/linalg/util/handle_functions.hpp"
#include "dca/parallel/util/get_workload.hpp"
#include "dca/phys/dca_step/cluster_solver/cluster_solver_id.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_accumulator.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_walker.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/qmci_autocorrelation_data.hpp"
#include "dca/phys/dca_step/cluster_solver/thread_task_handler.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/util/print_time.hpp"
#include "dca/distribution/dist_types.hpp"
#include "dca/math/function_transform/function_transform.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <class QmciSolver>
class StdThreadQmciClusterSolver : public QmciSolver {
public:
  using BaseClass = QmciSolver;
  using ThisType = StdThreadQmciClusterSolver<BaseClass>;
  static constexpr linalg::DeviceType device = QmciSolver::device;
  using Parameters = typename BaseClass::ParametersType;
  using Real = typename dca::config::McOptions::MC_REAL;
  using Scalar = typename dca::util::ScalarSelect<Real, Parameters::complex_g0>::type;
  using SignType = std::conditional_t<dca::util::IsComplex_t<Scalar>::value, Scalar, std::int8_t>;
  using Data = typename BaseClass::Data;
  using typename BaseClass::Concurrency;
  using typename BaseClass::Profiler;
  using typename BaseClass::Rng;
  using typename BaseClass::Accumulator;
  using Walker = stdthreadqmci::StdThreadQmciWalker<typename BaseClass::Walker, Data>;
  using SpGreensFunction = typename BaseClass::SpGreensFunction;
  using StdThreadAccumulatorType =
      stdthreadqmci::StdThreadQmciAccumulator<Accumulator, typename BaseClass::SpGreensFunction>;
  using MFunction = typename StdThreadAccumulatorType::MFunction;
  using MFunctionTime = typename StdThreadAccumulatorType::MFunctionTime;
  using MFunctionTimePair = typename StdThreadAccumulatorType::MFunctionTimePair;
  using FTauPair = typename StdThreadAccumulatorType::FTauPair;
  using PaddedTimeDmn = typename StdThreadAccumulatorType::PaddedTimeDmn;

protected:
  using BaseClass::accumulator_;

public:
  StdThreadQmciClusterSolver(Parameters& parameters_ref, Data& data_ref,
                             const std::shared_ptr<io::Writer<Concurrency>>& writer);

  void initialize(int dca_iteration);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

  void setSampleConfiguration(const io::Buffer& buff);

  struct MFuncAndSign {
    const MFunction& m_r_w;
    const SignType sign;
  };

  /** gets the MFunction which for CT-INT and CT-AUX is M_r_w
   *  we also need the sign since the MFunction is accumulated
   *  with as a product of fval and sign
   */
  MFuncAndSign getSingleMFunc(StdThreadAccumulatorType& accumulator) const {
    const MFunction& mfunc(accumulator.get_single_measurement_sign_times_MFunction());
    return {mfunc, accumulator.get_sign().getSign()};
  };

  struct MFuncTimeAndSign {
    const typename Accumulator::FTauPair& m_r_t;
    const SignType sign;
  };

  /** gets the MFunction in time domain  which for CT-INT and CT-AUX is M_r_t
   *  we also need the sign since the MFunction is accumulated
   *  with as a product of fval and sign
   */
  MFuncTimeAndSign getSingleMFuncTime(StdThreadAccumulatorType& accumulator) const {
    const FTauPair& mfunc(accumulator.get_single_measurement_sign_times_MFunction_time());
    return {mfunc, accumulator.get_sign().getSign()};
  };

  auto transformMFunction(const MFuncAndSign& mfs) const;
  auto computeSingleMeasurement_G_k_w(const SpGreensFunction& M_k_w) const;
  void logSingleMeasurement(StdThreadAccumulatorType& accumulator, int stamping_period,
                            bool log_MFunction, bool log_MFunctionTime) const;

private:
  void startWalker(int id);
  void startAccumulator(int id, const Parameters& parameters);
  void startWalkerAndAccumulator(int id, const Parameters& parameters);

  void initializeAndWarmUp(Walker& walker, int id, int walker_id);

  void readConfigurations();
  void writeConfigurations() const;

  void iterateOverLocalMeasurements(int walker_id, std::function<void(int, int, bool)>&& f);

  void printIntegrationMetadata() const;

  void finalizeWalker(Walker& walker, int walker_id);

private:
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

  dca::parallel::thread_traits::mutex_type mutex_merge_;
  dca::parallel::thread_traits::mutex_type mutex_queue_;
  dca::parallel::thread_traits::condition_variable_type queue_insertion_;

  std::vector<dca::io::Buffer> config_dump_;
  // stdthreadqmci::QmciAutocorrelationData<typename BaseClass::Walker> autocorrelation_data_;

  bool last_iteration_ = false;
  bool read_configuration_ = false;
  unsigned measurements_ = 0;
};

template <class QmciSolver>
StdThreadQmciClusterSolver<QmciSolver>::StdThreadQmciClusterSolver(
    Parameters& parameters_ref, Data& data_ref, const std::shared_ptr<io::Writer<Concurrency>>& writer)
    : BaseClass(parameters_ref, data_ref, writer),

      nr_walkers_(parameters_.get_walkers()),
      nr_accumulators_(parameters_.get_accumulators()),

      walker_fingerprints_(nr_walkers_, 0),
      accum_fingerprints_(nr_accumulators_, 0),

      thread_task_handler_(nr_walkers_, nr_accumulators_,
                           parameters_ref.shared_walk_and_accumulation_thread()),

      accumulators_queue_(),

      config_dump_(nr_walkers_)
      //autocorrelation_data_(parameters_, 0, BaseClass::g0_)
{
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
void StdThreadQmciClusterSolver<QmciSolver>::setSampleConfiguration(const io::Buffer& buff) {
  if (!parameters_.store_configuration())
    return;

  config_dump_.resize(nr_walkers_);
  for (auto& x : config_dump_)
    x = buff;

  read_configuration_ = true;
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

  std::vector<dca::parallel::thread_traits::future_type<void>> futures;

  dca::profiling::WallTime start_time;

  auto& pool = dca::parallel::ThreadPool::get_instance();
  for (int i = 0; i < thread_task_handler_.size(); ++i) {
    if (thread_task_handler_.getTask(i) == "walker")
      futures.emplace_back(pool.enqueue(&ThisType::startWalker, this, i));
    else if (thread_task_handler_.getTask(i) == "accumulator")
      futures.emplace_back(pool.enqueue(&ThisType::startAccumulator, this, i, parameters_));
    else if (thread_task_handler_.getTask(i) == "walker and accumulator")
      futures.emplace_back(pool.enqueue(&ThisType::startWalkerAndAccumulator, this, i, parameters_));
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

  if (parameters_.store_configuration()) {
    if (BaseClass::writer_) {
      if (BaseClass::writer_->isADIOS2()) {
        BaseClass::writer_->open_group("Configurations");
        BaseClass::writer_->execute("sample", config_dump_[0], true);
        BaseClass::writer_->close_group();
      }
      else if (concurrency_.id() == concurrency_.first()) {  // write one sample configuration.
        BaseClass::writer_->open_group("Configurations");
        BaseClass::writer_->execute("sample", config_dump_[0]);
        BaseClass::writer_->close_group();
      }
      read_configuration_ = true;
    }
  }
  QmciSolver::accumulator_.finalize();
  if (concurrency_.id() == concurrency_.first())
    std::cout << "accumulator finalized!" << std::endl;
}

template <class QmciSolver>
template <typename dca_info_struct_t>
double StdThreadQmciClusterSolver<QmciSolver>::finalize(dca_info_struct_t& dca_info_struct) {
  Profiler profiler(__FUNCTION__, "stdthread-MC-Integration", __LINE__);
  if (dca_iteration_ == parameters_.get_dca_iterations() - 1) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "Computing Error Bars." << std::endl;
    // For CTINT this is a no op.
    BaseClass::computeErrorBars();
  }

  // CTINT calculates its error here maybe
  double L2_Sigma_difference = QmciSolver::finalize(dca_info_struct);

  if (dca_iteration_ == parameters_.get_dca_iterations() - 1)
    writeConfigurations();

  // autocorrelation_data_.sumConcurrency(concurrency_);

  if (BaseClass::writer_ && *BaseClass::writer_ && concurrency_.id() == concurrency_.first()) {
    std::cout << "Writing actual run info\n";
    auto& writer = *BaseClass::writer_;
    writer.open_group("actual");
    int num_ranks = concurrency_.number_of_processors();
    writer.execute("ranks", num_ranks);
    std::vector<int> rank_measurements(num_ranks, 0);
    for (int ir = 0; ir < num_ranks; ++ir)
      rank_measurements[ir] = parallel::util::getWorkload(measurements_, num_ranks, ir);
    writer.execute("rank_measurements", rank_measurements);
    std::vector<int> thread_measurements(num_ranks * nr_walkers_, 0);
    for (int ir = 0; ir < num_ranks; ++ir)
      for (int iw = 0; iw < nr_walkers_; ++iw)
        thread_measurements[iw + ir * nr_walkers_] =
            parallel::util::getWorkload(rank_measurements[ir], nr_walkers_, iw);
    writer.execute("thread_measurements", thread_measurements);
    writer.close_group();

    // only CTAUX supports equal time accumulation.
    if constexpr (decltype(QmciSolver::accumulator_)::solver_id == ClusterSolverId::CT_AUX) {
      // This is a bit of a mess because we normally write the accumulator by writing the owning
      // integrator but currently this is expected to only happen 1 time per run.
      if (QmciSolver::accumulator_.perform_equal_time_accumulation()) {
        writer.open_group("CT-AUX-SOLVER-functions");
        QmciSolver::accumulator_.write(*BaseClass::writer_);
        writer.close_group();
      }
    }

    // Write and reset autocorrelation.
    std::cout << "Writing autocorrelation data\n";
    std::cout << "Autocorrelation incompatible with complex G0 and GPU";
    // autocorrelation_data_.write(*BaseClass::writer_, dca_iteration_);
  }
  // autocorrelation_data_.reset();

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

  auto walker_log = last_iteration_ ? BaseClass::writer_ : nullptr;
  Walker walker(parameters_, data_, rng_vector_[walker_index], concurrency_.get_id(), id,
                walker_log, BaseClass::g0_);

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
          // This can easily be out of sync with warmup ends
          if (print)
            walker.updateShell(meas_id, tot_meas);

          {
            Profiler profiler("stdthread-MC-walker waiting", "stdthread-MC-walker", __LINE__, id);
            acc_ptr = nullptr;

            // Wait for available accumulators.
            {
              dca::parallel::thread_traits::unique_lock lock(mutex_queue_);
              queue_insertion_.wait(lock, [&]() { return !accumulators_queue_.empty(); });
              acc_ptr = accumulators_queue_.front();
              accumulators_queue_.pop();
            }
          }
          acc_ptr->updateFrom(walker, concurrency_.id(), walker.get_thread_id(),
                              walker.get_meas_id(), last_iteration_);
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
    dca::parallel::thread_traits::scoped_lock lock(mutex_queue_);
    while (!accumulators_queue_.empty()) {
      accumulators_queue_.front()->notifyDone();
      accumulators_queue_.pop();
    }
  }

  finalizeWalker(walker, walker_index);

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
    throw std::runtime_error(
        "Non fix-meas-per-walker accumulation is suspect and disabled at this time.");
    // Perform the total number of loop with a shared atomic counter.
    // for (int meas_id = measurements_done_++; meas_id < n_local_meas; meas_id = measurements_done_++)
    //   f(meas_id, n_local_meas, print);
  }
}

template <class QmciSolver>
auto StdThreadQmciClusterSolver<QmciSolver>::transformMFunction(const MFuncAndSign& mfs) const {
  SpGreensFunction M;
  math::transform::FunctionTransform<typename QmciSolver::RDmn, typename QmciSolver::KDmn>::execute(
      mfs.m_r_w, M);
  M /= mfs.sign;
  return M;
}
// This repeated code is caused by the finalization of each single measurement M r w in the
// accumulators causing diversion from the normal algorithm. a way to merge this code should be
// found.
template <class QmciSolver>
auto StdThreadQmciClusterSolver<QmciSolver>::computeSingleMeasurement_G_k_w(
    const SpGreensFunction& M_k_w) const {
  SpGreensFunction G_k_w("G_k_w");
  QmciSolver::computeG_k_w(data_.G0_k_w_cluster_excluded, M_k_w, G_k_w);
  return G_k_w;
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::logSingleMeasurement(
    StdThreadAccumulatorType& accumulator_obj, int stamping_period, bool log_MFunction,
    bool log_MFunctionTime) const {
  if (accumulator_obj.get_meas_id() % stamping_period == 0) {
    if (log_MFunctionTime) {
      auto mfst = ThisType::getSingleMFuncTime(accumulator_obj);
      accumulator_obj.logPerConfigurationMFunctionTime(mfst.m_r_t, mfst.sign);
    }
    auto mfs = ThisType::getSingleMFunc(accumulator_obj);
    auto M_k_w = ThisType::transformMFunction(mfs);
    auto single_meas_G_k_w = ThisType::computeSingleMeasurement_G_k_w(M_k_w);
    if (log_MFunction)
      accumulator_obj.logPerConfigurationMFunction(M_k_w, mfs.sign);
    accumulator_obj.logPerConfigurationGreensFunction(single_meas_G_k_w);
    // We can remove this if we finally trust the accumulators to clear there single measurments.
    accumulator_obj.clearSingleMeasurement();
  }
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::startAccumulator(int id, const Parameters& parameters) {
  Profiler::start_threading(id);

  auto accumulator_log = last_iteration_ ? BaseClass::writer_ : nullptr;
  StdThreadAccumulatorType accumulator_obj(parameters_, data_, id, accumulator_log);

  accumulator_obj.initialize(dca_iteration_);

  std::unique_ptr<std::exception> exception_ptr;
  bool log_MFunction = parameters.per_measurement_MFunction();
  bool log_MFunctionTime = parameters.per_measurement_MFunction_time();
  auto stamping_period = parameters.stamping_period();
  try {
    while (true) {
      {
        dca::parallel::thread_traits::scoped_lock lock(mutex_queue_);
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
        if (accumulator_log && stamping_period)
          logSingleMeasurement(accumulator_obj, stamping_period, log_MFunction, log_MFunctionTime);
        accumulator_obj.finishMeasuring();
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

  assert(accumulator_obj.isMeasuring() == false);

  {
    dca::parallel::thread_traits::scoped_lock lock(mutex_merge_);
    accumulator_obj.sumTo(QmciSolver::accumulator_);

    // We want to do this here to avoid a race on writing the greens function.
    auto singleMeasurement_G0_write = last_iteration_ ? BaseClass::writer_ : nullptr;
  }

  accum_fingerprints_[thread_task_handler_.IDToAccumIndex(id)] = accumulator_obj.deviceFingerprint();
  Profiler::stop_threading(id);

  if (exception_ptr)
    throw(*exception_ptr);
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::startWalkerAndAccumulator(int id,
                                                                       const Parameters& parameters) {
  Profiler::start_threading(id);

  // Create and warm a walker.
  auto walker_log = BaseClass::writer_;
  Walker walker(parameters_, data_, rng_vector_[id], concurrency_.get_id(), id, walker_log,
                BaseClass::g0_);
  initializeAndWarmUp(walker, id, id);

  if (id == 0) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n\t\t warm-up ends\n" << std::endl;
  }

  auto accumulator_log = BaseClass::writer_;
  StdThreadAccumulatorType accumulator_obj(parameters_, data_, id, accumulator_log);
  accumulator_obj.initialize(dca_iteration_);

  std::unique_ptr<std::exception> current_exception;

  bool log_MFunction = parameters.per_measurement_MFunction();
  bool log_MFunctionTime = parameters.per_measurement_MFunction_time();
  auto stamping_period = parameters.stamping_period();

  try {
    iterateOverLocalMeasurements(id, [&](const int meas_id, const int n_meas, const bool print) {
      {
        Profiler profiler("Walker updating", "stdthread-MC", __LINE__, id);
        walker.doSweep();
      }
      {
        Profiler profiler("Accumulator measuring", "stdthread-MC", __LINE__, id);
        accumulator_obj.updateFrom(walker, concurrency_.get_id(), walker.get_thread_id(),
                                   walker.get_meas_id(), last_iteration_);
        accumulator_obj.measure();
        if (accumulator_log && stamping_period)
          logSingleMeasurement(accumulator_obj, stamping_period, log_MFunction, log_MFunctionTime);
        accumulator_obj.finishMeasuring();
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
  catch (...) {
    throw std::runtime_error("something mysterious  went wrong in walker thread!");
  }
			     
  ++walk_finished_;
  if (BaseClass::writer_ && BaseClass::writer_->isADIOS2())
    BaseClass::writer_->flush();

  {
    dca::parallel::thread_traits::scoped_lock lock(mutex_merge_);
    if (concurrency_.id() == concurrency_.first())
      std::cout << "Summing to accumulator --->";
    accumulator_obj.sumTo(QmciSolver::accumulator_);
  }
  if (concurrency_.id() == concurrency_.first())
    std::cout << " Done\n";

  finalizeWalker(walker, id);

  accum_fingerprints_[id] = accumulator_obj.deviceFingerprint();

  Profiler::stop_threading(id);

  if (current_exception)
    throw(*current_exception);
}

template <class QmciSolver>
void StdThreadQmciClusterSolver<QmciSolver>::finalizeWalker(Walker& walker, int walker_id) {
  config_dump_[walker_id] = walker.dumpConfig();
  walker_fingerprints_[walker_id] = walker.deviceFingerprint();
  // autocorrelation_data_ += walker;

  if (walker_id == 0 && concurrency_.id() == concurrency_.first()) {
    std::cout << "\n\t\t QMCI ends\n" << std::endl;
    walker.printSummary();
  }
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

    read_configuration_ = true;
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
