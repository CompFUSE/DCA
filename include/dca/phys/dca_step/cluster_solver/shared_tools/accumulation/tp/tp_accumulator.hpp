// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class computes the two particle Green's function from the walker's M matrix.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_HPP

#include <cmath>
#include <complex>
#include <memory>
#include <mutex>
#include <stdexcept>

#include "dca/config/config_defines.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrix_view.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/math/function_transform/special_transforms/space_transform_2D.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft_cpu.hpp"
#include "dca/phys/domains/cluster/momentum_exchange_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_exchange_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/models/traits.hpp"
#include "dca/phys/four_point_type.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

template <class Parameters, linalg::DeviceType device = linalg::CPU>
class TpAccumulator;

template <class Parameters>
class TpAccumulator<Parameters, linalg::CPU> {
public:
  using Real = typename Parameters::MC_measurement_scalar_type;

protected:
  using RDmn = typename Parameters::RClusterDmn;
  using KDmn = typename Parameters::KClusterDmn;
  using KExchangeDmn = func::dmn_0<domains::MomentumExchangeDomain>;
  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;
  using WDmn = func::dmn_0<domains::frequency_domain>;

public:
  using this_type = TpAccumulator<Parameters>;

  // Constructor:
  // In: G0: non interacting greens function.
  // In: pars: parameters object.
  // In: thread_id: thread id, only used by the profiler.
  TpAccumulator(
      const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
      const Parameters& pars, int thread_id = 0);

  // Resets the object between DCA iterations.
  void resetAccumulation(uint /*dca_loop*/ = 0);

  // Computes the two particles Greens function from the M matrix and accumulates it internally.
  // In: M_array: stores the M matrix for each spin sector.
  // In: configs: stores the walker's configuration for each spin sector.
  // In: sign: sign of the configuration.
  template <class Configuration, class Scalar>
  double accumulate(const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& M_pair,
                    const std::array<Configuration, 2>& configs, int sign);

  // Empty method for compatibility with GPU version.
  void finalize() {}

  // Returns the accumulated Green's function.
  const auto& get_sign_times_G4() const;

  // Sums the accumulated Green's function to the accumulated Green's function of other_acc.
  void sumTo(this_type& other_acc);

  void synchronizeCopy() {}

protected:
  using Profiler = typename Parameters::profiler_type;

  using WTpDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using WTpPosDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT_POSITIVE>>;
  using WTpExtDmn = func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED>>;
  using WTpExtPosDmn = func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED_POSITIVE>>;
  using WExchangeDmn = func::dmn_0<domains::FrequencyExchangeDomain>;

  using Complex = std::complex<Real>;
  using SpGreenFunction =
      func::function<Complex, func::dmn_variadic<BDmn, BDmn, SDmn, KDmn, KDmn, WTpExtPosDmn, WTpExtDmn>>;

  using TpDomain =
      func::dmn_variadic<BDmn, BDmn, BDmn, BDmn, KDmn, KDmn, KExchangeDmn, WTpDmn, WTpDmn, WExchangeDmn>;
  using TpGreenFunction = func::function<Complex, TpDomain>;
  using Matrix = linalg::Matrix<Complex, linalg::CPU>;

  void initializeG0();

  double computeG();

  void computeGMultiband(int s, int k1, int k2, int w1, int w2);

  void computeGSingleband(int s, int k1, int k2, int w1, int w2);

  void getGMultiband(int s, int k1, int k2, int w1, int w2, Matrix& G, Complex beta = 0) const;

  Complex getGSingleband(int s, int k1, int k2, int w1, int w2) const;

  template <class Configuration, typename Scalar>
  double computeM(const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& M_pair,
                  const std::array<Configuration, 2>& configs);

  double updateG4();

  void inline updateG4Atomic(Complex* G4_ptr, int s_a, int k1_a, int k2_a, int w1_a, int w2_a,
                             int s_b, int k1_b, int k2_b, int w1_b, int w2_b, Real alpha,
                             bool cross_legs);

  void inline updateG4SpinDifference(Complex* G4_ptr, int sign, int k1_a, int k2_a, int w1_a,
                                     int w2_a, int k1_b, int k2_b, int w1_b, int w2_b, Real alpha,
                                     bool cross_legs);

protected:
  const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>* const G0_ptr_ =
      nullptr;

  const int thread_id_;

  const FourPointType mode_ = NONE;
  const Real beta_ = -1;
  constexpr static int n_bands_ = Parameters::model_type::BANDS;

  constexpr static bool non_density_density_ =
      models::has_non_density_interaction<typename Parameters::lattice_type>::value;
  CachedNdft<Real, RDmn, WTpExtDmn, WTpExtPosDmn, linalg::CPU, non_density_density_> ndft_obj_;

  SpGreenFunction G_;

  std::unique_ptr<TpGreenFunction> G4_;

  func::function<Complex, func::dmn_variadic<BDmn, BDmn, SDmn, KDmn, WTpExtDmn>> G0_;

  int sign_;

  const int extension_index_offset_ = -1;
  const int n_pos_frqs_ = -1;

private:
  // work spaces for computeGMultiband.
  Matrix G0_M_, G_a_, G_b_;
};

template <class Parameters>
TpAccumulator<Parameters, linalg::CPU>::TpAccumulator(
    const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
    const Parameters& pars, const int thread_id)
    : G0_ptr_(&G0),
      thread_id_(thread_id),
      mode_(pars.get_four_point_type()),
      beta_(pars.get_beta()),
      extension_index_offset_((WTpExtDmn::dmn_size() - WTpDmn::dmn_size()) / 2),
      n_pos_frqs_(WTpExtPosDmn::dmn_size()),
      G0_M_(n_bands_),
      G_a_(n_bands_),
      G_b_(n_bands_) {
  if (WDmn::dmn_size() < WTpExtDmn::dmn_size())
    throw(std::logic_error("The number of single particle frequencies is too small."));
  initializeG0();
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::CPU>::resetAccumulation(uint /*dca_loop*/) {
  G4_.reset(new TpGreenFunction("G4"));
  initializeG0();
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::CPU>::initializeG0() {
  const int sp_index_offset = (WDmn::dmn_size() - WTpExtDmn::dmn_size()) / 2;

  for (int w = 0; w < WTpExtDmn::dmn_size(); ++w) {
    assert(std::abs(WTpExtDmn::get_elements()[w] - WDmn::get_elements()[w + sp_index_offset]) < 1e-3);
    for (int k = 0; k < KDmn::dmn_size(); ++k)
      for (int s = 0; s < 2; ++s)
        for (int b2 = 0; b2 < n_bands_; ++b2)
          for (int b1 = 0; b1 < n_bands_; ++b1)
            G0_(b1, b2, s, k, w) = (*G0_ptr_)(b1, s, b2, s, k, w + sp_index_offset);
  }
}

template <class Parameters>
template <class Configuration, class Scalar>
double TpAccumulator<Parameters, linalg::CPU>::accumulate(
    const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& M_pair,
    const std::array<Configuration, 2>& configs, const int sign) {
  Profiler profiler("accumulate", "tp-accumulation", __LINE__, thread_id_);
  double gflops(0.);
  if (!(configs[0].size() + configs[1].size()))  // empty config
    return gflops;

  sign_ = sign;
  gflops += computeM(M_pair, configs);
  gflops += computeG();
  gflops += updateG4();
  return gflops;
}

template <class Parameters>
template <class Configuration, class Scalar>
double TpAccumulator<Parameters, linalg::CPU>::computeM(
    const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& M_pair,
    const std::array<Configuration, 2>& configs) {
  double gflops(0.);

  func::function<Complex, func::dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, WTpExtPosDmn, WTpExtDmn>> M_r_r_w_w;

  for (int spin = 0; spin < SDmn::dmn_size(); ++spin) {
    Profiler prf_a("Frequency FT", "tp-accumulation", __LINE__, thread_id_);
    if (not configs[spin].size())
      continue;
    gflops += ndft_obj_.execute(configs[spin], M_pair[spin], M_r_r_w_w, spin);
  }

  Profiler prf_b("Space FT", "tp-accumulation", __LINE__, thread_id_);
  // TODO: add the gflops here.
  math::transform::SpaceTransform2D<RDmn, KDmn, Real>::execute(M_r_r_w_w, G_);

  return gflops;
}

template <class Parameters>
double TpAccumulator<Parameters, linalg::CPU>::computeG() {
  Profiler prf("ComputeG", "tp-accumulation", __LINE__, thread_id_);
  for (int w2 = 0; w2 < WTpExtDmn::dmn_size(); ++w2)
    for (int w1 = 0; w1 < WTpExtPosDmn::dmn_size(); ++w1)
      for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
        for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1)
          for (int s = 0; s < 2; ++s)
            switch (n_bands_) {
              case 1:
                computeGSingleband(s, k1, k2, w1, w2);
                break;
              default:
                computeGMultiband(s, k1, k2, w1, w2);
            }
  //  INTERNAL: the additional flops for w1==w2 are ignored.
  const Real flops = 8 * std::pow(n_bands_, 3) * WTpExtPosDmn::dmn_size() * WTpExtDmn::dmn_size() *
                     std::pow(KDmn::dmn_size(), 2) * 2;
  return 1e-9 * flops;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::CPU>::computeGSingleband(const int s, const int k1,
                                                                const int k2, const int w1,
                                                                const int w2) {
  assert(w1 < WTpExtPosDmn::dmn_size());
  assert(w2 < WTpExtDmn::dmn_size());

  const Complex G0_w1 = G0_(0, 0, s, k1, w1 + n_pos_frqs_);
  const Complex G0_w2 = G0_(0, 0, s, k2, w2);
  const Complex M_val = G_(0, 0, s, k1, k2, w1, w2);

  if (k2 == k1 && w2 == w1 + n_pos_frqs_)
    G_(0, 0, s, k1, k2, w1, w2) = -G0_w1 * M_val * G0_w2 + G0_w1 * beta_;
  else
    G_(0, 0, s, k1, k2, w1, w2) = -G0_w1 * M_val * G0_w2;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::CPU>::computeGMultiband(const int s, const int k1,
                                                               const int k2, const int w1,
                                                               const int w2) {
  assert(w1 < WTpExtPosDmn::dmn_size());
  assert(w2 < WTpExtDmn::dmn_size());

  const linalg::MatrixView<Complex, linalg::CPU> G0_w1(&G0_(0, 0, s, k1, w1 + n_pos_frqs_),
                                                       n_bands_, n_bands_);
  const linalg::MatrixView<Complex, linalg::CPU> G0_w2(&G0_(0, 0, s, k2, w2), n_bands_, n_bands_);
  linalg::MatrixView<Complex, linalg::CPU> M_matrix(&G_(0, 0, s, k1, k2, w1, w2), n_bands_);

  // G(w1, w2) <- -G0(w1) M(w1, w2) G0(w2)
  linalg::matrixop::gemm(G0_w1, M_matrix, G0_M_);
  linalg::matrixop::gemm(Complex(-1), G0_M_, G0_w2, Complex(0), M_matrix);

  // G(w1, w2) += \delta(w1, w2) \delta(k1,k2) G0(w1)
  if (G0_w1.ptr() == G0_w2.ptr()) {
    for (int b2 = 0; b2 < n_bands_; ++b2)
      for (int b1 = 0; b1 < n_bands_; ++b1)
        M_matrix(b1, b2) += G0_w1(b1, b2) * beta_;
  }
}

template <class Parameters>
std::complex<typename TpAccumulator<Parameters, linalg::CPU>::Real> TpAccumulator<
    Parameters, linalg::CPU>::getGSingleband(const int s, const int k1, const int k2, const int w1,
                                             const int w2) const {
  const int w2_ext = w2 + extension_index_offset_;
  const int w1_ext = w1 + extension_index_offset_;
  auto minus_w1 = [=](const int w) { return n_pos_frqs_ - 1 - w; };
  auto minus_w2 = [=](const int w) { return 2 * n_pos_frqs_ - 1 - w; };
  auto plus_w1 = [=](const int w) { return w - n_pos_frqs_; };

  if (w1_ext >= n_pos_frqs_)
    return G_(0, 0, s, k1, k2, plus_w1(w1_ext), w2_ext);
  else
    return std::conj(G_(0, 0, s, k1, k2, minus_w1(w1_ext), minus_w2(w2_ext)));
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::CPU>::getGMultiband(int s, int k1, int k2, int w1, int w2,
                                                           Matrix& G, const Complex beta) const {
  const int w2_ext = w2 + extension_index_offset_;
  const int w1_ext = w1 + extension_index_offset_;
  auto minus_w1 = [=](const int w) { return n_pos_frqs_ - 1 - w; };
  auto minus_w2 = [=](const int w) { return 2 * n_pos_frqs_ - 1 - w; };
  auto plus_w1 = [=](const int w) { return w - n_pos_frqs_; };

  if (w1_ext >= n_pos_frqs_) {
    const Complex* const G_ptr = &G_(0, 0, s, k1, k2, plus_w1(w1_ext), w2_ext);
    for (int b2 = 0; b2 < n_bands_; ++b2)
      for (int b1 = 0; b1 < n_bands_; ++b1)
        G(b1, b2) = beta * G(b1, b2) + G_ptr[b1 + b2 * n_bands_];
  }
  else {
    const Complex* const G_ptr = &G_(0, 0, s, k1, k2, minus_w1(w1_ext), minus_w2(w2_ext));
    for (int b2 = 0; b2 < n_bands_; ++b2)
      for (int b1 = 0; b1 < n_bands_; ++b1)
        G(b1, b2) = beta * G(b1, b2) + std::conj(G_ptr[b1 + b2 * n_bands_]);
  }
}

template <class Parameters>
double TpAccumulator<Parameters, linalg::CPU>::updateG4() {
  // G4 is stored with the following band convention:
  // b1 ------------------------ b3
  //        |           |
  //        |           |
  //        |           |
  // b2 ------------------------ b4
  Profiler profiler("updateG4", "tp-accumulation", __LINE__, thread_id_);
  if (!G4_)
    G4_.reset(new TpGreenFunction("G4"));
  double flops(0);

  auto momentum_sum = [](const int k, const int q) { return KDmn::parameter_type::add(k, q); };
  auto q_minus_k = [](const int k, const int q) { return KDmn::parameter_type::subtract(k, q); };
  // Returns the index of the exchange frequency w_ex plus the Matsubara frequency with index w.
  auto w_plus_w_ex = [](const int w, const int w_ex) { return w + w_ex; };
  // Returns the index of the exchange frequency w_ex minus the Matsubara frequency with index w.
  auto w_ex_minus_w = [](const int w, const int w_ex) { return w_ex + WTpDmn::dmn_size() - 1 - w; };

  const Real sign_over_2 = 0.5 * sign_;

  const double flops_update_atomic = 3 * std::pow(n_bands_, 4);
  const double flops_update_spin_diff = flops_update_atomic + 2 * std::pow(n_bands_, 2);
  const int n_loops = WTpDmn::dmn_size() * WTpDmn::dmn_size() * KDmn::dmn_size() * KDmn::dmn_size();

  const auto& exchange_frq = domains::FrequencyExchangeDomain::get_elements();
  const auto& exchange_mom = domains::MomentumExchangeDomain::get_elements();

  switch (mode_) {
    case PARTICLE_HOLE_MAGNETIC:
      //       G4(k1, k2, k_ex) = 1/2  (s1 * s2) <c^+(k1 + k_ex, s1) c(k1, s1)
      //                       c^+(k2, s2) c(k2 + k_ex, s2)>
      //                     = 1/2 (s1 * s2) <G(k1, k1 + k_ex, s1) G(k2 + k_ex, k2, s2) - (s1 ==s2)
      //                       G(k2 + k_ex, k1 + k_ex, s1) G(k1, k2, s1)>
      for (int w_ex_idx = 0; w_ex_idx < exchange_frq.size(); ++w_ex_idx) {
        const int w_ex = exchange_frq[w_ex_idx];
        for (int w2 = 0; w2 < WTpDmn::dmn_size(); ++w2)
          for (int w1 = 0; w1 < WTpDmn::dmn_size(); ++w1)
            for (int k_ex_idx = 0; k_ex_idx < exchange_mom.size(); ++k_ex_idx) {
              const int k_ex = exchange_mom[k_ex_idx];
              for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
                for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
                  Complex* const G4_ptr = &(*G4_)(0, 0, 0, 0, k1, k2, k_ex_idx, w1, w2, w_ex_idx);
                  updateG4SpinDifference(G4_ptr, -1, k1, momentum_sum(k1, k_ex), w1,
                                         w_plus_w_ex(w1, w_ex), momentum_sum(k2, k_ex), k2,
                                         w_plus_w_ex(w2, w_ex), w2, sign_over_2, false);
                  for (int s = 0; s < 2; ++s)
                    updateG4Atomic(G4_ptr, s, k1, k2, w1, w2, s, momentum_sum(k2, k_ex),
                                   momentum_sum(k1, k_ex), w_plus_w_ex(w2, w_ex),
                                   w_plus_w_ex(w1, w_ex), -sign_over_2, true);
                }
            }
      }
      flops += n_loops * (flops_update_spin_diff + 2 * flops_update_atomic);
      break;

    case PARTICLE_HOLE_CHARGE:
      // G4(k1, k2, k_ex) += 1/2  <c^+(k1 + k_ex, s1) c(k1, s1) c^+(k2, s2) c(k2 + k_ex, s2)> =
      //                = 1/2 <G(k1, k1 + k_ex, s1) G(k2 + k_ex, k2, s2) -
      //                  (s1 ==s2) G(k2 + k_ex, k1 + k_ex, s1) G(k1, k2, s1)>
      for (int w_ex_idx = 0; w_ex_idx < exchange_frq.size(); ++w_ex_idx) {
        const int w_ex = exchange_frq[w_ex_idx];
        for (int w2 = 0; w2 < WTpDmn::dmn_size(); ++w2)
          for (int w1 = 0; w1 < WTpDmn::dmn_size(); ++w1)
            for (int k_ex_idx = 0; k_ex_idx < exchange_mom.size(); ++k_ex_idx) {
              const int k_ex = exchange_mom[k_ex_idx];
              for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
                for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
                  Complex* const G4_ptr = &(*G4_)(0, 0, 0, 0, k1, k2, k_ex_idx, w1, w2, w_ex_idx);
                  updateG4SpinDifference(G4_ptr, 1, k1, momentum_sum(k1, k_ex), w1,
                                         w_plus_w_ex(w1, w_ex), momentum_sum(k2, k_ex), k2,
                                         w_plus_w_ex(w2, w_ex), w2, sign_over_2, false);
                  for (int s = 0; s < 2; ++s)
                    updateG4Atomic(G4_ptr, s, k1, k2, w1, w2, s, momentum_sum(k2, k_ex),
                                   momentum_sum(k1, k_ex), w_plus_w_ex(w2, w_ex),
                                   w_plus_w_ex(w1, w_ex), -sign_over_2, true);
                }
            }
      }

      flops += n_loops * (flops_update_spin_diff + 2 * flops_update_atomic);
      break;

    case PARTICLE_HOLE_TRANSVERSE:
      // G4 = 1/2 \sum_s <c^+(k1 + k_ex, s) c(k1, -s) c^+(k2, -s) c(k2 + k_ex, s)>
      //    = -1/2 \sum_s G( k2 + k_ex, k1 + k_ex, s) G4(k1, k2, -s)
      for (int w_ex_idx = 0; w_ex_idx < exchange_frq.size(); ++w_ex_idx) {
        const int w_ex = exchange_frq[w_ex_idx];
        for (int w2 = 0; w2 < WTpDmn::dmn_size(); ++w2)
          for (int w1 = 0; w1 < WTpDmn::dmn_size(); ++w1)
            for (int k_ex_idx = 0; k_ex_idx < exchange_mom.size(); ++k_ex_idx) {
              const int k_ex = exchange_mom[k_ex_idx];
              for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
                for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
                  Complex* const G4_ptr = &(*G4_)(0, 0, 0, 0, k1, k2, k_ex_idx, w1, w2, w_ex_idx);
                  for (int s = 0; s < 2; ++s)
                    updateG4Atomic(G4_ptr, s, k1, k2, w1, w2, not s, momentum_sum(k2, k_ex),
                                   momentum_sum(k1, k_ex), w_plus_w_ex(w2, w_ex),
                                   w_plus_w_ex(w1, w_ex), -sign_over_2, true);
                }
            }
      }

      flops += n_loops * 2 * flops_update_atomic;
      break;

    case PARTICLE_PARTICLE_UP_DOWN:
      // G4 = 1/2 \sum_s <c^+(k_ex-k1, s) c^+(k1, -s) c(k2, -s) c(k_ex-k2, s)>
      //    = 1/2 \sum_s G(k_ex-k2, k_ex-k1, s) G4(k2, k1, -s)
      for (int w_ex_idx = 0; w_ex_idx < exchange_frq.size(); ++w_ex_idx) {
        const int w_ex = exchange_frq[w_ex_idx];
        for (int w2 = 0; w2 < WTpDmn::dmn_size(); ++w2)
          for (int w1 = 0; w1 < WTpDmn::dmn_size(); ++w1)
            for (int k_ex_idx = 0; k_ex_idx < exchange_mom.size(); ++k_ex_idx) {
              const int k_ex = exchange_mom[k_ex_idx];
              for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
                for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
                  Complex* const G4_ptr = &(*G4_)(0, 0, 0, 0, k1, k2, k_ex_idx, w1, w2, w_ex_idx);
                  for (int s = 0; s < 2; ++s)
                    updateG4Atomic(G4_ptr, s, k1, k2, w1, w2, !s, q_minus_k(k1, k_ex),
                                   q_minus_k(k2, k_ex), w_ex_minus_w(w1, w_ex),
                                   w_ex_minus_w(w2, w_ex), sign_over_2, false);
                }
            }
      }

      flops += n_loops * 2 * flops_update_atomic;
      break;
    default:
      throw(std::logic_error("Non supported tp mode."));
  }

  return 1e-9 * flops;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::CPU>::updateG4Atomic(
    Complex* G4_ptr, const int s_a, const int k1_a, const int k2_a, const int w1_a, const int w2_a,
    const int s_b, const int k1_b, const int k2_b, const int w1_b, const int w2_b, const Real alpha,
    const bool cross_legs) {
  // This function performs the following update for each band:
  //
  // G4(k1, k2, w1, w2) += alpha * G(s_a, k1_a, k2_a, w1_a, w2_a) * G(s_b, k1_b, k2_b, w1_b, w2_b)
  // INTERNAL: would use __restrict__ pointer make sense?
  if (n_bands_ == 1) {
    *G4_ptr += alpha * getGSingleband(s_a, k1_a, k2_a, w1_a, w2_a) *
               getGSingleband(s_b, k1_b, k2_b, w1_b, w2_b);
  }
  else {
    getGMultiband(s_a, k1_a, k2_a, w1_a, w2_a, G_a_);
    getGMultiband(s_b, k1_b, k2_b, w1_b, w2_b, G_b_);

    if (!cross_legs)
      for (int b4 = 0; b4 < n_bands_; ++b4)
        for (int b3 = 0; b3 < n_bands_; ++b3)
          for (int b2 = 0; b2 < n_bands_; ++b2)
            for (int b1 = 0; b1 < n_bands_; ++b1) {
              *G4_ptr += alpha * G_a_(b1, b3) * G_b_(b2, b4);
              ++G4_ptr;
            }
    else
      for (int b4 = 0; b4 < n_bands_; ++b4)
        for (int b3 = 0; b3 < n_bands_; ++b3)
          for (int b2 = 0; b2 < n_bands_; ++b2)
            for (int b1 = 0; b1 < n_bands_; ++b1) {
              *G4_ptr += alpha * G_a_(b1, b4) * G_b_(b2, b3);
              ++G4_ptr;
            }
  }
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::CPU>::updateG4SpinDifference(
    Complex* G4_ptr, const int sign, const int k1_a, const int k2_a, const int w1_a, const int w2_a,
    const int k1_b, const int k2_b, const int w1_b, const int w2_b, const Real alpha,
    const bool cross_legs) {
  // This function performs the following update for each band:
  //
  // G4(k1, k2, w1, w2) += alpha * (G(up, k1_a, k2_a, w1_a, w2_a) * G(up, k1_b, k2_b, w1_b, w2_b)
  //                       + sign * G(down, k1_a, k2_a, w1_a, w2_a) * G(down, k1_b, k2_b, w1_b,
  //                       w2_b)
  if (n_bands_ == 1) {
    *G4_ptr += alpha * (getGSingleband(0, k1_a, k2_a, w1_a, w2_a) +
                        Complex(sign) * getGSingleband(1, k1_a, k2_a, w1_a, w2_a)) *
               (getGSingleband(0, k1_b, k2_b, w1_b, w2_b) +
                Complex(sign) * getGSingleband(1, k1_b, k2_b, w1_b, w2_b));
  }
  else {
    getGMultiband(0, k1_a, k2_a, w1_a, w2_a, G_a_);
    getGMultiband(1, k1_a, k2_a, w1_a, w2_a, G_a_, sign);
    getGMultiband(0, k1_b, k2_b, w1_b, w2_b, G_b_);
    getGMultiband(1, k1_b, k2_b, w1_b, w2_b, G_b_, sign);

    if (!cross_legs)
      for (int b4 = 0; b4 < n_bands_; ++b4)
        for (int b3 = 0; b3 < n_bands_; ++b3)
          for (int b2 = 0; b2 < n_bands_; ++b2)
            for (int b1 = 0; b1 < n_bands_; ++b1) {
              *G4_ptr += alpha * G_a_(b1, b3) * G_b_(b2, b4);
              ++G4_ptr;
            }
    else
      for (int b4 = 0; b4 < n_bands_; ++b4)
        for (int b3 = 0; b3 < n_bands_; ++b3)
          for (int b2 = 0; b2 < n_bands_; ++b2)
            for (int b1 = 0; b1 < n_bands_; ++b1) {
              *G4_ptr += alpha * G_a_(b1, b4) * G_b_(b2, b3);
              ++G4_ptr;
            }
  }
}

template <class Parameters>
const auto& TpAccumulator<Parameters, linalg::CPU>::get_sign_times_G4() const {
  if (!G4_)
    throw(std::logic_error("There is no G4 stored in this class."));
  return *G4_;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::CPU>::sumTo(this_type& other_one) {
  if (!G4_)
    throw(std::logic_error("There is no G4 stored in this class."));

  if (!other_one.G4_)
    other_one.G4_.reset(new TpGreenFunction("G4"));
  *(other_one.G4_) += *G4_;

  G4_.release();
}

}  // accumulator
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_HPP
