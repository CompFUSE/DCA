// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class measures the correlation of G(r = 0, t = 0).

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_ACCUMULATOR_TP_TIME_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_ACCUMULATOR_TP_TIME_ACCUMULATOR_HPP

#include <mutex>

#include "dca/math/statistics/autocorrelation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_matrix_configuration.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/linalg/matrix_view.hpp"
#include "dca/linalg/multi_vector.hpp"

#ifdef DCA_HAVE_CUDA
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration_manager.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/kernels_interface.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation_gpu.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

using dca::linalg::CPU;
using dca::linalg::GPU;
using dca::linalg::DeviceType;

template <class Parameters, typename Real, DeviceType device>
class TimeCorrelator {
public:
  using Profiler = typename Parameters::profiler_type;
  using Concurrency = typename Parameters::concurrency_type;

  using TDmn = func::dmn_0<domains::time_domain>;
  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using RDmn = typename Parameters::RClusterDmn;

  template <DeviceType conf_device>
  using Config = linalg::MultiVector<conf_device, Real, int, int>;

  TimeCorrelator(const Parameters& parameters_, int id);

  template <class WalkerConfig, typename RealInp>
  void compute_G_r_t(const std::array<dca::linalg::Matrix<RealInp, device>, 2>& M,
                     const std::array<WalkerConfig, 2>& configs, int sign);

  static void setG0(const ctint::G0Interpolation<device, Real>& g0) {
    g0_ = &g0;
  }

  auto& getCorrelators() {
    return autocorrelations_;
  }

  auto& operator+=(const TimeCorrelator& rhs);
  void sumConcurrency(const Concurrency& concurrency);

  void reset();

  float get_FLOPs() const;

private:
  void initializeFixedConfiguration();

  template <class WalkerConfig>
  void uploadConfig(const std::array<WalkerConfig, 2>& config, int s);

  void computeG0(linalg::Matrix<Real, GPU>& G0_mat, const Config<GPU>& config_l,
                 const Config<GPU>& config_r);

  //  template <typename = std::enable_if_t<device == CPU>>
  void computeG0(linalg::Matrix<Real, CPU>& G0_mat, const Config<CPU>& config_l,
                 const Config<CPU>& config_r);

  const Parameters& parameters_;
  const Concurrency& concurrency_;

  int thread_id_;
  float flops_ = 0;
  constexpr static int n_bands_ = Parameters::bands;

  constexpr static int n_correlators_ =
      n_bands_ * (n_bands_ + 1) / 2;  // Number of independent entries in a band-band matrix.

  static const inline ctint::G0Interpolation<device, Real>* g0_ = nullptr;

  const linalg::util::CudaStream& stream_;

  linalg::Matrix<Real, device> G0_l_;
  linalg::Matrix<Real, device> G0_r_;
  linalg::Matrix<Real, device> G0_M_;
  std::array<linalg::Matrix<Real, device>, 2> G_;

  std::array<linalg::Matrix<Real, CPU>, 2> G_host_;

  static inline std::unique_ptr<Config<device>> fixed_config_;
  std::array<Config<device>, 2> m_r_dev_config_;
  std::array<Config<device>, 2> m_l_dev_config_;
  std::array<Config<CPU>, 2> m_l_host_config_;
  std::array<Config<CPU>, 2> m_r_host_config_;

  std::array<math::statistics::Autocorrelation<Real>, n_correlators_> autocorrelations_;
};

template <class Parameters, typename Real, DeviceType device>
TimeCorrelator<Parameters, Real, device>::TimeCorrelator(const Parameters& parameters_ref, int id)
    : parameters_(parameters_ref),
      concurrency_(parameters_.get_concurrency()),

      thread_id_(id),
      stream_(linalg::util::getStream(thread_id_, 0)) {
  static std::once_flag flag;
  std::call_once(flag, [&]() {
    initializeFixedConfiguration();
#ifdef DCA_HAVE_CUDA
    // TODO: rename helper to ClusterHelper.
    ctint::CtintHelper::set<RDmn, BDmn>();
#endif
  });

  for (auto& correlator : autocorrelations_) {
    correlator.resize(parameters_.get_time_correlation_window());
  }
}

template <class Parameters, typename Real, DeviceType device>
void TimeCorrelator<Parameters, Real, device>::initializeFixedConfiguration() {
  fixed_config_ = std::make_unique<Config<device>>(n_bands_);

  Config<CPU> host_config(n_bands_);
  Real* tau = host_config.template get<0>();
  int* b = host_config.template get<1>();
  int* r = host_config.template get<2>();

  for (int i = 0; i < n_bands_; ++i) {
    tau[i] = 0;
    r[i] = 0;
    b[i] = i;
  }

  fixed_config_->setAsync(host_config, stream_);
  stream_.sync();
}

template <class Parameters, typename Real, DeviceType device>
template <class WalkerConfig, typename RealInp>
void TimeCorrelator<Parameters, Real, device>::compute_G_r_t(
    const std::array<dca::linalg::Matrix<RealInp, device>, 2>& M,
    const std::array<WalkerConfig, 2>& configs, int sign) {
  if (!g0_) {
    throw(std::runtime_error("G0 is not set."));
  }

  // Upload state
  constexpr int n_electron_spins = 1;  // Compute only on up-up sector.

  for (int s = 0; s < n_electron_spins; ++s)
    uploadConfig(configs, s);

  for (int s = 0; s < n_electron_spins; ++s) {
    G_[s].resizeNoCopy(fixed_config_->size());

    if (M[s].nrRows() == 0) {
      G_host_[s].setToZero(stream_);
      continue;
    }

    G0_M_.resizeNoCopy(std::make_pair(fixed_config_->size(), m_l_dev_config_[s].size()));

    computeG0(G0_l_, *fixed_config_, m_l_dev_config_[s]);
    flops_ += 9. * G0_l_.nrRows() * G0_l_.nrRows();
    linalg::matrixop::gemm(G0_l_, M[s], G0_M_, thread_id_, 0);
    flops_ += 2. * G0_l_.nrRows() * G0_l_.nrCols() * G0_M_.nrCols();

    computeG0(G0_r_, m_r_dev_config_[s], *fixed_config_);
    flops_ += 9. * G0_r_.nrRows() * G0_r_.nrRows();
    linalg::matrixop::gemm(G0_M_, G0_r_, G_[s], thread_id_, 0);
    flops_ += 2. * G0_M_.nrRows() * G0_M_.nrCols() * G_[s].nrCols();

    G_host_[s].setAsync(G_[s], stream_);
  }

  linalg::util::syncStream(thread_id_, 0);

  // Note: the constant addendum G0(fixed_configuration) is skipped.

  int linindex = 0;
  for (int b1 = 0; b1 < n_bands_; ++b1)
    for (int b2 = b1; b2 < n_bands_; ++b2, ++linindex) {
      autocorrelations_[linindex].addSample(sign * G_host_[0](b1, b2));
    }
}

template <class Parameters, typename Real, DeviceType device>
template <class WalkerConfig>
void TimeCorrelator<Parameters, Real, device>::uploadConfig(const std::array<WalkerConfig, 2>& configs,
                                                            const int s) {
  const auto& sector = configs[s];
  m_l_host_config_[s].resizeNoCopy(sector.size());
  m_r_host_config_[s].resizeNoCopy(sector.size());

  Real* t_l = m_l_host_config_[s].template get<0>();
  int* b_l = m_l_host_config_[s].template get<1>();
  int* r_l = m_l_host_config_[s].template get<2>();

  Real* t_r = m_r_host_config_[s].template get<0>();
  int* b_r = m_r_host_config_[s].template get<1>();
  int* r_r = m_r_host_config_[s].template get<2>();

  for (int i = 0; i < sector.size(); ++i) {
    t_l[i] = t_r[i] = sector[i].get_tau();
    // Note: M is an inverse G0, hence left and right are swapped.
    b_l[i] = sector[i].get_right_band();
    r_l[i] = sector[i].get_right_site();
    b_r[i] = sector[i].get_left_band();
    r_r[i] = sector[i].get_left_site();
  }

  m_l_dev_config_[s].setAsync(m_l_host_config_[s], stream_);
  m_r_dev_config_[s].setAsync(m_r_host_config_[s], stream_);
}

template <class Parameters, typename Real, DeviceType device>
float TimeCorrelator<Parameters, Real, device>::get_FLOPs() const {
  return flops_;
}

template <class Parameters, typename Real, DeviceType device>
auto& TimeCorrelator<Parameters, Real, device>::operator+=(
    const TimeCorrelator<Parameters, Real, device>& rhs) {
  for (int i = 0; i < autocorrelations_.size(); ++i) {
    autocorrelations_[i] += rhs.autocorrelations_[i];
  }

  return *this;
}

template <class Parameters, typename Real, DeviceType device>
void TimeCorrelator<Parameters, Real, device>::sumConcurrency(const Concurrency& concurrency) {
  // TODO: delay.
  for (int i = 0; i < autocorrelations_.size(); ++i) {
    autocorrelations_[i].sumConcurrency(concurrency);
  }
}

template <class Parameters, typename Real, DeviceType device>
void TimeCorrelator<Parameters, Real, device>::reset() {
  for (auto& c : autocorrelations_)
    c.reset();
}

#ifdef DCA_HAVE_CUDA
template <class Parameters, typename Real, DeviceType device>
void TimeCorrelator<Parameters, Real, device>::computeG0(linalg::Matrix<Real, GPU>& G0_mat,
                                                         const Config<GPU>& config_l,
                                                         const Config<GPU>& config_r) {
  G0_mat.resizeNoCopy(std::make_pair(config_l.size(), config_r.size()));

  linalg::MatrixView<Real, device> G0_view(G0_mat);
  details::computeG0(G0_view, *g0_, config_l.template get<0>(), config_l.template get<1>(),
                     config_l.template get<2>(), config_r.template get<0>(),
                     config_r.template get<1>(), config_r.template get<2>(), stream_);
}
#endif

template <class Parameters, typename Real, DeviceType device>
void TimeCorrelator<Parameters, Real, device>::computeG0(linalg::Matrix<Real, CPU>& G0_mat,
                                                         const Config<CPU>& config_l,
                                                         const Config<CPU>& config_r) {
  G0_mat.resizeNoCopy(std::make_pair(config_l.size(), config_r.size()));

  using BDmn = func::dmn_0<domains::electron_band_domain>;
  func::dmn_variadic<BDmn, BDmn, RDmn> labels;

  const Real* t_l = config_l.template get<0>();
  const int* b_l = config_l.template get<1>();
  const int* r_l = config_l.template get<2>();

  const Real* t_r = config_r.template get<0>();
  const int* b_r = config_r.template get<1>();
  const int* r_r = config_r.template get<2>();

  for (int j = 0; j < G0_mat.nrCols(); ++j)
    for (int i = 0; i < G0_mat.nrRows(); ++i) {
      const Real tau = t_l[i] - t_r[j];
      const int r = RDmn::parameter_type::subtract(r_r[j], r_l[i]);
      const int label = labels(b_l[i], b_r[j], r);

      G0_mat(i, j) = (*g0_)(tau, label);
    }
}

}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_ACCUMULATOR_TP_TIME_ACCUMULATOR_HPP
