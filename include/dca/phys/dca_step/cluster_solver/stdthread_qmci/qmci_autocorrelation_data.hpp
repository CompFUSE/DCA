// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// A std::thread jacket that measures correlations in the MC walker independent of the MC method.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_QMCI_AUTOCORRELATION_DATA_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_QMCI_AUTOCORRELATION_DATA_HPP

#include <mutex>

#include "dca/io/writer.hpp"
#include "dca/math/statistics/autocorrelation.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/time_correlator.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_walker.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace stdthreadqmci {
// dca::phys::solver::stdthreadqmci::

/** provides autocorrelation data functionality to walkers as base class
 *  unfortunately does other things as well and causes unclear code involving walkers.
 *  \todo rethink this poor design.
 */
template <class Walker>
class QmciAutocorrelationData {
  using Parameters = typename Walker::parameters_type;
  using Data = DcaData<Parameters>;
  using Concurrency = typename Parameters::concurrency_type;
  using Real = typename dca::config::McOptions::MC_REAL;
  using Scalar = typename dca::util::ScalarSelect<Real,Parameters::complex_g0>::type;

  constexpr static auto device = Walker::device;
  constexpr static int bands = Parameters::bands;

public:
  QmciAutocorrelationData(const Parameters& parameters, int thread_id,
                          G0Interpolation<device, Scalar>& g0);
  virtual ~QmciAutocorrelationData() = default;

  void accumulateAutocorrelation(Walker& walker);

  void sumConcurrency(const Concurrency& concurrency);

  // Accumulate the data from the other object. This method is thread safe.
  QmciAutocorrelationData& operator+=(const QmciAutocorrelationData& other);

  void write(dca::io::Writer<Concurrency>& writer, int dca_loop);

  void reset();

  void markThermalized();

private:
  const unsigned autocorrelation_window_;
  const bool accumulate_G_;

  std::array<dca::linalg::Matrix<Scalar, device>, 2> m_correlator_;

  TimeCorrelator<Parameters, Scalar, device> time_correlator_;
  math::statistics::Autocorrelation<int> order_correlator_;
  math::statistics::Autocorrelation<Scalar> weight_correlator_;

  // Store MC weights for each chain
  using SignType = dca::util::SignType<Scalar>;
  std::vector<std::vector<SignType>> signs_;
  std::vector<std::vector<double>> weights_;
  std::vector<std::vector<unsigned long>> steps_;
  std::vector<unsigned> thermalization_step_;
};

template <class Walker>
QmciAutocorrelationData<Walker>::QmciAutocorrelationData(const Parameters& parameters,
                                                         const int thread_id,
                                                         G0Interpolation<device, Scalar>& g0)
    : autocorrelation_window_(parameters.get_time_correlation_window()),
      accumulate_G_(parameters.compute_G_correlation()),
      time_correlator_(parameters, thread_id, g0),
      order_correlator_(autocorrelation_window_),
      weight_correlator_(autocorrelation_window_) {}

template <class Walker>
void QmciAutocorrelationData<Walker>::write(io::Writer<Concurrency>& writer, int dca_loop) {
  // Write MC weights
  writer.open_group("Configurations");
  writer.open_group("MC-weight-samples");
  writer.open_group("iteration " + std::to_string(dca_loop));

  writer.execute("steps", steps_);
  writer.execute("signs", signs_);
  writer.execute("weights", weights_);
  writer.execute("thermalization-step", thermalization_step_);

  writer.close_group();
  writer.close_group();
  writer.close_group();

  if (!autocorrelation_window_)
    return;

  writer.open_group("Autocorrelation");
  writer.open_group("iteration " + std::to_string(dca_loop));

  // Write G(t = 0).
  if (accumulate_G_) {
    writer.open_group("G_t0");

    linalg::Matrix<Scalar, linalg::CPU> g_corr(bands, "autocorr");
    linalg::Matrix<Scalar, linalg::CPU> g_stdev(bands, "stdev");
    linalg::Matrix<Scalar, linalg::CPU> g_mean(bands, "mean");

    int lindex = 0;
    for (int b1 = 0; b1 < bands; ++b1)
      for (int b2 = b1; b2 < bands; ++b2, ++lindex) {
        auto& correlator = time_correlator_.getCorrelators()[lindex];
        g_corr(b1, b2) = g_corr(b2, b1) = correlator.computeAutocorrelationTime();
        g_stdev(b1, b2) = g_stdev(b2, b1) = correlator.getStdev();
        g_mean(b1, b2) = g_mean(b2, b1) = correlator.getMean();
      }

    writer.execute(g_corr);
    writer.execute(g_stdev);
    writer.execute(g_mean);
    writer.close_group();
  }

  // Write expansion order and weight
  auto write_correlator = [&](auto& correlator, const std::string& name) {
    writer.open_group(name);
    writer.execute("autocorrelation", correlator.computeAutocorrelationTime());
    writer.execute("stdev", correlator.getStdev());
    writer.execute("mean", correlator.getMean());
    writer.close_group();
  };

  write_correlator(order_correlator_, "expansion order");
  write_correlator(weight_correlator_, "MC weight");

  writer.close_group();
  writer.close_group();
}

template <class Walker>
void QmciAutocorrelationData<Walker>::accumulateAutocorrelation(Walker& walker) {
  if (weights_.size() == 0) {
    weights_.resize(1);
    signs_.resize(1);
    steps_.resize(1);
  }

  weights_.back().push_back(walker.get_MC_log_weight());
  signs_.back().push_back(walker.get_sign());
  steps_.back().push_back(walker.get_steps());

  if (autocorrelation_window_ && walker.is_thermalized()) {
    if (accumulate_G_) {
      walker.computeM(m_correlator_);
      time_correlator_.compute_G_r_t(m_correlator_, walker.get_matrix_configuration(),
                                     walker.get_sign());
    }

    order_correlator_.addSample(walker.get_configuration().size());

    weight_correlator_.addSample(walker.get_MC_log_weight());
  }
}

template <class Walker>
void QmciAutocorrelationData<Walker>::markThermalized() {
  assert(thermalization_step_.size() == 0);
  thermalization_step_.push_back(steps_.back().back());
}

template <class Walker>
QmciAutocorrelationData<Walker>& QmciAutocorrelationData<Walker>::operator+=(
    const QmciAutocorrelationData<Walker>& other) {
  static std::mutex mutex;
  std::unique_lock<std::mutex> lock(mutex);

  // Collect weight measurements.
  signs_.insert(signs_.end(), other.signs_.begin(), other.signs_.end());
  weights_.insert(weights_.end(), other.weights_.begin(), other.weights_.end());
  steps_.insert(steps_.end(), other.steps_.begin(), other.steps_.end());
  thermalization_step_.insert(thermalization_step_.end(), other.thermalization_step_.begin(),
                              other.thermalization_step_.end());

  if (autocorrelation_window_) {
    if (accumulate_G_)
      time_correlator_ += other.time_correlator_;
    order_correlator_ += other.order_correlator_;
    weight_correlator_ += other.weight_correlator_;
  }

  return *this;
}

template <class Walker>
void QmciAutocorrelationData<Walker>::sumConcurrency(const Concurrency& concurrency) {
  if (autocorrelation_window_) {
    if (accumulate_G_)
      time_correlator_.sumConcurrency(concurrency);
    order_correlator_.sumConcurrency(concurrency);
    weight_correlator_.sumConcurrency(concurrency);
  }

  // Don't communicate MC weights. (Too much data).
}

template <class Walker>
void QmciAutocorrelationData<Walker>::reset() {
  weights_.clear();
  signs_.clear();
  steps_.clear();
  thermalization_step_.clear();

  time_correlator_.reset();
  order_correlator_.reset();
  weight_correlator_.reset();
}

// No additional ops specialization for SS-HYB
template <linalg::DeviceType device, class Parameters, class Data>
class QmciAutocorrelationData<cthyb::SsCtHybWalker<device, Parameters, Data>> {
  using Concurrency = typename Parameters::concurrency_type;
  using Walker = cthyb::SsCtHybWalker<device, Parameters, Data>;
public:
  using Real = typename Walker::Scalar;

  QmciAutocorrelationData(const Parameters&, int, G0Interpolation<device, Real>& g0) : g0_(g0) {}
  void sumConcurrency(const Concurrency&) {}
  QmciAutocorrelationData& operator+=(const QmciAutocorrelationData&) {
    return *this;
  }
  static void write([[maybe_unused]] io::Writer<Concurrency>& writer, int) {}
  void reset() {}
protected:
  G0Interpolation<device, Real>& g0_;
};

}  // namespace stdthreadqmci
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_QMCI_AUTOCORRELATION_DATA_HPP
