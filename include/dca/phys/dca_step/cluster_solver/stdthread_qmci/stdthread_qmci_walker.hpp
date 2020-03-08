// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// A std::thread jacket that measures correlations in the MC walker independent of the MC method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_WALKER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_WALKER_HPP

#include <mutex>

#include "dca/math/statistics/autocorrelation.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/time_correlator.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_walker.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace stdthreadqmci {
// dca::phys::solver::stdthreadqmci::

template <class QmciWalker>
class StdThreadQmciWalker : public QmciWalker {
  using ThisType = StdThreadQmciWalker<QmciWalker>;
  using Parameters = typename QmciWalker::parameters_type;
  using Data = DcaData<Parameters>;
  using Concurrency = typename Parameters::concurrency_type;
  using Rng = typename Parameters::random_number_generator;
  using Real = typename QmciWalker::Scalar;

  constexpr static auto device = QmciWalker::device;
  constexpr static int bands = Parameters::bands;

public:
  StdThreadQmciWalker(const Parameters& parameters_ref, Data& data_ref, Rng& rng, int id,
                      io::HDF5Writer* writer);
  ~StdThreadQmciWalker();

  void doSweep();

  static void write(io::HDF5Writer& writer);
  static void sumConcurrency(const Concurrency& concurrency);

private:
  const unsigned autocorrelation_window_;
  const unsigned stamping_period_;
  int thread_id_;

  std::array<dca::linalg::Matrix<Real, device>, 2> m_correlator_;

  TimeCorrelator<Parameters, Real, device> time_correlator_;
  math::statistics::Autocorrelation<int> order_correlator_;

  io::HDF5Writer* const writer_ = nullptr;
  std::size_t meas_id_ = 0;

  static inline std::unique_ptr<TimeCorrelator<Parameters, Real, device>> common_time_correlator_;
  static inline math::statistics::Autocorrelation<int> common_order_correlator_;
};

template <class QmciWalker>
StdThreadQmciWalker<QmciWalker>::StdThreadQmciWalker(const Parameters& parameters, Data& data_ref,
                                                     Rng& rng, const int id, io::HDF5Writer* writer)
    : QmciWalker(parameters, data_ref, rng, id),
      autocorrelation_window_(parameters.get_time_correlation_window()),
      stamping_period_(parameters.stamping_period()),
      thread_id_(id),
      time_correlator_(parameters, thread_id_),
      order_correlator_(autocorrelation_window_),
      writer_(writer) {
  static std::once_flag flag;
  std::call_once(flag, [&]() {
    if (autocorrelation_window_) {
      common_time_correlator_ =
          std::make_unique<TimeCorrelator<Parameters, Real, device>>(parameters, id);
      common_order_correlator_.resize(autocorrelation_window_);
    }
  });
}

template <class QmciWalker>
void StdThreadQmciWalker<QmciWalker>::doSweep() {
  QmciWalker::doSweep();

  if (autocorrelation_window_ && QmciWalker::is_thermalized()) {
    QmciWalker::computeM(m_correlator_);
    time_correlator_.compute_G_r_t(m_correlator_, QmciWalker::get_matrix_configuration(),
                                   QmciWalker::get_sign());
    order_correlator_.addSample(QmciWalker::get_configuration().size());
  }

  ++meas_id_;
  const bool print_to_log = writer_ && static_cast<bool>(*writer_);  // File exists and it is open.
  if (print_to_log && stamping_period_ && (meas_id_ % stamping_period_) == 0) {
    const std::string stamp_name =
        "w_" + std::to_string(thread_id_) + "_step_" + std::to_string(meas_id_);

    writer_->lock();

    writer_->open_group("Configurations");
    auto& config = QmciWalker::get_configuration();
    config.write(*writer_, stamp_name);
    writer_->close_group();

    writer_->unlock();
  }
}

template <class QmciWalker>
StdThreadQmciWalker<QmciWalker>::~StdThreadQmciWalker() {
  if (autocorrelation_window_) {
    static std::mutex mutex;
    std::unique_lock<std::mutex> lock(mutex);

    time_correlator_.sumTo(*common_time_correlator_);
    order_correlator_.sumTo(common_order_correlator_);
  }
}

template <class QmciWalker>
void StdThreadQmciWalker<QmciWalker>::sumConcurrency(const Concurrency& concurrency) {
  common_time_correlator_->sumConcurrency(concurrency);
  common_order_correlator_.sumConcurrency(concurrency);
}

template <class QmciWalker>
void StdThreadQmciWalker<QmciWalker>::write(io::HDF5Writer& writer) {
  linalg::Matrix<double, linalg::CPU> g_corr(bands, "G_t0_autocorr");
  linalg::Matrix<double, linalg::CPU> g_stdev(bands, "G_t0_stdev");

  int lindex = 0;
  for (int b1 = 0; b1 < bands; ++b1)
    for (int b2 = b1; b2 < bands; ++b2, ++lindex) {
      auto& correlator = common_time_correlator_->getCorrelators()[lindex];
      g_corr(b1, b2) = g_corr(b2, b1) = correlator.computeAutocorrelationTime();
      g_stdev(b1, b2) = g_stdev(b2, b1) = correlator.getStdev();
    }

  writer.execute(g_corr);
  writer.execute(g_stdev);
  common_time_correlator_->reset();

  writer.execute("order_autocorrelation", common_order_correlator_.computeAutocorrelationTime());
  writer.execute("order_stdev", common_order_correlator_.getStdev());
  common_order_correlator_.reset();
}

// No additional ops specialization for SS-HYB
template <linalg::DeviceType device, class Parameters, class Data>
class StdThreadQmciWalker<cthyb::SsCtHybWalker<device, Parameters, Data>>
    : public cthyb::SsCtHybWalker<device, Parameters, Data> {
  using QmciWalker = cthyb::SsCtHybWalker<device, Parameters, Data>;
  using ThisType = StdThreadQmciWalker<QmciWalker>;
  using Rng = typename Parameters::random_number_generator;

public:
  StdThreadQmciWalker(const Parameters& parameters_ref, Data& data_ref, Rng& rng, int id,
                      io::HDF5Writer* /*writer*/)
      : QmciWalker(parameters_ref, data_ref, rng, id) {}

  static void write(io::HDF5Writer& writer) {}
};

}  // namespace stdthreadqmci
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_STDTHREAD_QMCI_WALKER_HPP
