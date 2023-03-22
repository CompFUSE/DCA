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
public:
  using ThisType = StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>;
  using Parameters = typename QmciAccumulator::ParametersType;
  using Real = typename dca::config::McOptions::MC_REAL;
  using Scalar = typename dca::util::ScalarSelect<Real,Parameters::complex_g0>::type;
  using SignType = std::conditional_t<dca::util::IsComplex_t<Scalar>::value, Scalar, std::int8_t>;
  using Concurrency = typename Parameters::concurrency_type;
  using Data = typename QmciAccumulator::DataType;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using CDA = ClusterDomainAliases<Parameters::lattice_type::DIMENSION>;
  using KDmn = typename CDA::KClusterDmn;
  using RDmn = typename CDA::RClusterDmn;
  using MFunction = typename QmciAccumulator::MFunction;
  using MFunctionTime = typename QmciAccumulator::MFunctionTime;
  using MFunctionTimePair = typename QmciAccumulator::MFunctionTimePair;
  using FTauPair = typename QmciAccumulator::FTauPair;
  using PaddedTimeDmn = typename QmciAccumulator::PaddedTimeDmn;
  
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

  void logPerConfigurationMFunction(const SpGreensFunction&, const SignType sign) const;
  void logPerConfigurationMFunctionTime(const typename QmciAccumulator::FTauPair&,
                                        const SignType sign) const;

  void measure();

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sumTo(QmciAccumulator& other);

  void clearSingleMeasurementM_r_t() {
    QmciAccumulator::clearSingleMeasurement();
  }
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

  std::size_t get_meas_id() const {
    return meas_id_;
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
  // Temp hack to deal with vexing ss_ct_hyb build issue
  MFunction dummy_mfunc;
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
void StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>::logPerConfigurationMFunction(
    const SpGreensFunction& mfunc, const SignType sign) const {
  const bool print_to_log = writer_ && static_cast<bool>(*writer_);  // File exists and it is open. \todo possibly this should always be true
  if (print_to_log && stamping_period_ && (meas_id_ % stamping_period_) == 0) {
    if (writer_ && (writer_->isADIOS2() || concurrency_id_ == 0)) {
      const std::string stamp_name = "r_" + std::to_string(concurrency_id_) + "_meas_" +
                                     std::to_string(meas_id_) + "_w_" +
                                     std::to_string(walker_thread_id_);
      auto signFreeMFunc = mfunc;
      signFreeMFunc /= -sign;
      writer_->lock();
      writer_->open_group("STQW_Configurations");
      writer_->open_group(stamp_name);
      writer_->execute("MFunction", signFreeMFunc);
      writer_->execute("sign", sign);
      writer_->close_group();
      writer_->close_group();
      writer_->unlock();
    }
  }
}

template <class QmciAccumulator, class SpGreensFunction>
void StdThreadQmciAccumulator<QmciAccumulator, SpGreensFunction>::logPerConfigurationMFunctionTime(
    const typename QmciAccumulator::FTauPair& mfunc, const SignType sign) const {
  const bool print_to_log = writer_ && static_cast<bool>(*writer_);  // File exists and it is open.
  if (print_to_log && stamping_period_ && (meas_id_ % stamping_period_) == 0) {
    if (writer_ && (writer_->isADIOS2() || concurrency_id_ == 0)) {
      const std::string stamp_name = "r_" + std::to_string(concurrency_id_) + "_meas_" +
                                     std::to_string(meas_id_) + "_w_" +
                                     std::to_string(walker_thread_id_);

      using MFTauSpin =
	func::function<typename QmciAccumulator::Scalar, func::dmn_variadic<SDmn, PaddedTimeDmn>>;

      MFTauSpin mfunc_func;

      for (int i_spin = 0; i_spin < 2; ++i_spin) {
        std::copy_n(mfunc[0].data(), mfunc[0].size(), mfunc_func.values() + mfunc_func.subind_2_linind(i_spin,0));
      }
      mfunc_func /= -sign;
      writer_->lock();
      writer_->open_group("STQW_Configurations");
      writer_->open_group(stamp_name);
      writer_->execute("MFunctionTime", mfunc_func);
      writer_->execute("sign", sign);
      writer_->close_group();
      writer_->close_group();
      writer_->unlock();
    }
  }
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
