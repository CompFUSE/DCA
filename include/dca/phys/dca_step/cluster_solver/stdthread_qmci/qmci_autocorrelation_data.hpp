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

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_QMCI_AUTOCORRELATION_DATA_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_QMCI_AUTOCORRELATION_DATA_HPP

#include <mutex>

#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/math/statistics/autocorrelation.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/time_correlator.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_walker.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace stdthreadqmci {
// dca::phys::solver::stdthreadqmci::

template <class Walker>
class QmciAutocorrelationData {
  using Parameters = typename Walker::parameters_type;
  using Data = DcaData<Parameters>;
  using Concurrency = typename Parameters::concurrency_type;
  using Real = typename Walker::Scalar;

  constexpr static auto device = Walker::device;
  constexpr static int bands = Parameters::bands;

public:
  QmciAutocorrelationData(const Parameters& parameters, int thread_id);
  virtual ~QmciAutocorrelationData() = default;

  void accumulateAutocorrelation(Walker& walker);

  void sumConcurrency(const Concurrency& concurrency);

  // Accumulate the data from the other object. This method is thread safe.
  QmciAutocorrelationData& operator+=(const QmciAutocorrelationData& other);

  void write(dca::io::HDF5Writer& writer, int dca_loop);

  void reset();

private:
  const unsigned autocorrelation_window_;
  const double log_beta_;
  const bool accumulate_G_;

  std::array<dca::linalg::Matrix<Real, device>, 2> m_correlator_;

  TimeCorrelator<Parameters, Real, device> time_correlator_;
  math::statistics::Autocorrelation<int> order_correlator_;
  math::statistics::Autocorrelation<Real> energy_correlator_;
};

template <class Walker>
QmciAutocorrelationData<Walker>::QmciAutocorrelationData(const Parameters& parameters,
                                                         const int thread_id)
    : autocorrelation_window_(parameters.get_time_correlation_window()),
      log_beta_(std::log(parameters.get_beta())),
      accumulate_G_(parameters.compute_G_correlation()),
      time_correlator_(parameters, thread_id),
      order_correlator_(autocorrelation_window_),
      energy_correlator_(autocorrelation_window_) {}

template <class Walker>
void QmciAutocorrelationData<Walker>::write(io::HDF5Writer& writer, int dca_loop) {
  if (!autocorrelation_window_)
    return;

  writer.open_group("Autocorrelation");
  writer.open_group("iteration " + std::to_string(dca_loop));

  // Write G(t = 0).
  if (accumulate_G_) {
    writer.open_group("G_t0");

    linalg::Matrix<Real, linalg::CPU> g_corr(bands, "autocorr");
    linalg::Matrix<Real, linalg::CPU> g_stdev(bands, "stdev");
    linalg::Matrix<Real, linalg::CPU> g_mean(bands, "mean");

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

  // Write expansion order and energy
  auto write_correlator = [&](auto& correlator, const std::string& name) {
    writer.open_group(name);
    writer.execute("autocorrelation", correlator.computeAutocorrelationTime());
    writer.execute("stdev", correlator.getStdev());
    writer.execute("mean", correlator.getMean());
    writer.close_group();
  };

  write_correlator(order_correlator_, "expansion order");
  write_correlator(energy_correlator_, "energy");

  writer.close_group();
  writer.close_group();
}

template <class Walker>
void QmciAutocorrelationData<Walker>::accumulateAutocorrelation(Walker& walker) {
  if (autocorrelation_window_ && walker.is_thermalized()) {
    if (accumulate_G_) {
      walker.computeM(m_correlator_);
      time_correlator_.compute_G_r_t(m_correlator_, walker.get_matrix_configuration(),
                                     walker.get_sign());
    }

    order_correlator_.addSample(walker.get_configuration().size());

    const Real energy = -(walker.get_MC_log_weight() - log_beta_);
    energy_correlator_.addSample(energy);
  }
}

template <class Walker>
QmciAutocorrelationData<Walker>& QmciAutocorrelationData<Walker>::operator+=(
    const QmciAutocorrelationData<Walker>& other) {
  if (autocorrelation_window_) {
    static std::mutex mutex;
    std::unique_lock<std::mutex> lock(mutex);

    if (accumulate_G_)
      time_correlator_ += other.time_correlator_;
    order_correlator_ += other.order_correlator_;
    energy_correlator_ += other.energy_correlator_;
  }

  return *this;
}

template <class Walker>
void QmciAutocorrelationData<Walker>::sumConcurrency(const Concurrency& concurrency) {
  if (autocorrelation_window_) {
    if (accumulate_G_)
      time_correlator_.sumConcurrency(concurrency);
    order_correlator_.sumConcurrency(concurrency);
    energy_correlator_.sumConcurrency(concurrency);
  }
}

template <class Walker>
void QmciAutocorrelationData<Walker>::reset() {
  time_correlator_.reset();
  order_correlator_.reset();
  energy_correlator_.reset();
}

// No additional ops specialization for SS-HYB
template <linalg::DeviceType device, class Parameters, class Data>
class QmciAutocorrelationData<cthyb::SsCtHybWalker<device, Parameters, Data>> {
  using Concurrency = typename Parameters::concurrency_type;

public:
  QmciAutocorrelationData(const Parameters&, int) {}
  void sumConcurrency(const Concurrency&) {}
  QmciAutocorrelationData& operator+=(const QmciAutocorrelationData&) {
    return *this;
  }
  static void write(io::HDF5Writer&, int) {}
  void reset() {}
};

}  // namespace stdthreadqmci
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_STDTHREAD_QMCI_QMCI_AUTOCORRELATION_DATA_HPP
