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
// A std::thread jacket that implements a MC accumulator independent of the MC method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP

#include <atomic>
#include <queue>
#include <stdexcept>

#include "dca/io/writer.hpp"
#include "dca/config/threading.hpp"
#include "dca/math/function_transform/function_transform.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace stdthreadqmci {
// dca::phys::solver::stdthreadqmci::

template <class QmciAccumulator, class SpGreensFunction>
class StdThreadQmciAccumulator : public QmciAccumulator {
  using ThisType = StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>;
  using Parameters = typename QmciAccumulator::ParametersType;
  using Data = typename QmciAccumulator::DataType;
  using CDA = ClusterDomainAliases<Parameters::lattice_type::DIMENSION>;
  using KDmn = typename CDA::KClusterDmn;
  using RDmn = typename CDA::RClusterDmn;

public:
  StdThreadQmciAccumulator(const Parameters& parameters_ref, Data& data_ref, int id,
                           std::shared_ptr<io::Writer<Concurrency>> writer);

  ~StdThreadQmciAccumulator();

  using QmciAccumulator::finalize;
  using QmciAccumulator::initialize;

  template <typename Walker>
  void updateFrom(Walker& walker, int concurrency_id, int walker_thread_id, std::size_t meas_id,
                  bool last_iteration);

  void waitForQmciWalker();

  void logPerConfigurationGreensFunction(const SpGreensFunction&) const;

  void measure();

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sumTo(QmciAccumulator& other);

  const auto& get_single_measurment_sign_times_M_r_w() {
    return QmciAccumulator::get_single_measurment_sign_times_M_r_w();
  };
  // Signals that this object will not need to perform any more accumulation.
  void notifyDone();

  bool done() const {
    return done_;
  }

  bool isMeasuring() const {
    return measuring_;
  }

  void finishMeasuring() {
    measuring_ = false;
  }

private:
  int thread_id_;
  bool measuring_;
  int concurrency_id_;
  int walker_thread_id_;
  std::size_t meas_id_;
  bool last_iteration_;
  std::atomic<bool> done_;
  dca::parallel::thread_traits::condition_variable_type start_measuring_;
  dca::parallel::thread_traits::mutex_type mutex_accumulator_;
  const unsigned stamping_period_;
  std::shared_ptr<io::Writer<Concurrency>> writer_;
  const Data& data_ref_;
};

template <class QmciAccumulator, class SpGreensFunction>
StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>::StdThreadQmciAccumulator(
    const Parameters& parameters_ref, Data& data_ref, const int id,
    std::shared_ptr<io::Writer<Concurrency>> writer)
    : QmciAccumulator(parameters_ref, data_ref, id),
      thread_id_(id),
      measuring_(false),
      done_(false),
      stamping_period_(parameters_ref.stamping_period()),
      writer_(writer),
      data_ref_(data_ref) {}

template <class QmciAccumulator, class SpGreensFunction>
StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>::~StdThreadQmciAccumulator() {}

template <class QmciAccumulator, class SpGreensFunction>
template <typename Walker>
void StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>::updateFrom(Walker& walker,
                                                                             int concurrency_id,
                                                                             int walker_thread_id,
                                                                             std::size_t meas_id,
                                                                             bool last_iteration) {
  {
    // take a lock and keep it until it goes out of scope
    dca::parallel::thread_traits::unique_lock lock(mutex_accumulator_);
    if (measuring_)
      throw std::logic_error(__FUNCTION__);

    QmciAccumulator::updateFrom(walker);
    measuring_ = true;
    concurrency_id_ = concurrency_id;
    walker_thread_id_ = walker_thread_id;
    meas_id_ = meas_id;
    last_iteration_ = last_iteration;
  }

  start_measuring_.notify_one();
}

template <class QmciAccumulator, class SpGreensFunction>
void StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>::waitForQmciWalker() {
  dca::parallel::thread_traits::unique_lock lock(mutex_accumulator_);
  start_measuring_.wait(lock, [this]() { return measuring_ || done_; });
}

template <class QmciAccumulator, class SpGreensFunction>
void StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>::measure() {
  dca::parallel::thread_traits::scoped_lock lock(mutex_accumulator_);

  if (done_)
    return;
  assert(measuring_);

  QmciAccumulator::measure();
}

template <class QmciAccumulator, class SpGreensFunction>
void StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>::logPerConfigurationGreensFunction(
    const SpGreensFunction& spf) const {
  const bool print_to_log = writer_ && static_cast<bool>(*writer_);  // File exists and it is open.
  if (print_to_log && stamping_period_ && (meas_id_ % stamping_period_) == 0) {
    if (writer_ && (writer_->isADIOS2() || concurrency_id_ == 0)) {
      const std::string stamp_name = "r_" + std::to_string(concurrency_id_) + "_meas_" +
                                     std::to_string(meas_id_) + "_w_" +
                                     std::to_string(walker_thread_id_);
      writer_->lock();
      writer_->open_group("STQW_Configurations");
      writer_->open_group(stamp_name);
      writer_->execute("G_k_w", spf);
      writer_->close_group();
      writer_->close_group();
      writer_->unlock();
    }
  }
}

template <class QmciAccumulator, class SpGreensFunction>
void StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>::sumTo(QmciAccumulator& other) {
  dca::parallel::thread_traits::scoped_lock lock(mutex_accumulator_);
  QmciAccumulator::sumTo(other);
}

template <class QmciAccumulator, class SpGreensFunction>
void StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>::notifyDone() {
  done_ = true;
  start_measuring_.notify_one();
}

}  // namespace stdthreadqmci
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_ACCUMULATOR_HPP
