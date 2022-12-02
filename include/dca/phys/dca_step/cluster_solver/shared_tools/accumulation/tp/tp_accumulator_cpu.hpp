// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//

#ifndef DCA_TP_ACCUMULATOR_CPU_HPP
#define DCA_TP_ACCUMULATOR_CPU_HPP

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "dca/config/config_defines.hpp"
#include "dca/distribution/dist_types.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrix_view.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/linalg/util/gpu_stream.hpp"
#include "dca/math/function_transform/special_transforms/space_transform_2D.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft_cpu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_base.hpp"

#include "dca/phys/domains/cluster/momentum_exchange_domain.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_base.hpp"
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

template <class Parameters, DistType DT>
class TpAccumulator<Parameters, DT, linalg::CPU> : public TpAccumulatorBase<Parameters, DT> {
public:
  using Base = TpAccumulatorBase<Parameters, DT>;

  using typename Base::NuDmn;
  using typename Base::Real;
  using typename Base::Complex;
  using typename Base::RDmn;
  using typename Base::KDmn;
  using typename Base::WDmn;
  using typename Base::WTpDmn;
  using typename Base::WTpExtDmn;
  using typename Base::WTpExtPosDmn;
  using typename Base::WTpPosDmn;
  using typename Base::BDmn;
  using typename Base::SDmn;
  using typename Base::TpGreensFunction;
  
protected:
  using Base::non_density_density_;
  using Base::n_bands_;
  using Base::extension_index_offset_;
  using Base::n_pos_frqs_;
  using Base::G4_;
  using Base::channels_;
  using Base::G0_;
  using Base::G0_ptr_;
  using Base::G_;
  using Base::beta_;
  
  using Profiler = typename Parameters::profiler_type;
  using Base::thread_id_;
  
  using Matrix = linalg::Matrix<Complex, linalg::CPU>;

public:
  // Constructor:
  // In: G0: non interacting greens function.
  // In: pars: parameters object.
  // In: thread_id: thread id, only used by the profiler.
  TpAccumulator(
      const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
      const Parameters& pars, int thread_id = 0);

  // Resets the object between DCA iterations.
  void resetAccumulation(unsigned int /*dca_loop*/ = 0);

  // Computes the two particles Greens function from the M matrix and accumulates it internally.
  // In: M_array: stores the M matrix for each spin sector.
  // In: configs: stores the walker's configuration for each spin sector.
  // In: sign: sign of the configuration.
  template <class Configuration, typename RealIn>
  double accumulate(const std::array<linalg::Matrix<RealIn, linalg::CPU>, 2>& M_pair,
                    const std::array<Configuration, 2>& configs, int sign);

  // Empty method for compatibility with GPU version.
  void finalize() {}

  // Empty method for compatibility with GPU version.
  void ringG() {}

  // Returns the accumulated Green's function.
  const std::vector<TpGreensFunction>& get_G4() const;

  // FOR TESTING: Returns the accumulated Green's function.
  std::vector<TpGreensFunction>& get_nonconst_G4();

  // Sums the accumulated Green's function to the accumulated Green's function of other_acc.
  void sumTo(TpAccumulator& other_acc);

  void synchronizeCopy() {}

  template <class T>
  void syncStreams(const T&) {}

  std::size_t deviceFingerprint() const {
    return 0;
  }
  static std::size_t staticDeviceFingerprint() {
    return 0;
  }

  const linalg::util::GpuStream* get_stream() const {
    static const dca::linalg::util::GpuStream mock_stream;
    return &mock_stream;
  }

protected:
  double computeG();

  void computeGMultiband(int s, int k1, int k2, int w1, int w2);

  void computeGSingleband(int s, int k1, int k2, int w1, int w2);

  void getGMultiband(int s, int k1, int k2, int w1, int w2, Matrix& G, Complex beta = 0) const;

  auto getGSingleband(int s, int k1, int k2, int w1, int w2) -> Complex const;

  template <class Configuration, typename RealIn>
  float computeM(const std::array<linalg::Matrix<RealIn, linalg::CPU>, 2>& M_pair,
                 const std::array<Configuration, 2>& configs);

  double updateG4(int channel_id);

  void inline updateG4Atomic(Complex* G4_ptr, int s_a, int k1_a, int k2_a, int w1_a, int w2_a,
                             int s_b, int k1_b, int k2_b, int w1_b, int w2_b, Real alpha,
                             bool cross_legs);

  void inline updateG4SpinDifference(Complex* G4_ptr, int sign, int k1_a, int k2_a, int w1_a,
                                     int w2_a, int k1_b, int k2_b, int w1_b, int w2_b, Real alpha,
                                     bool cross_legs);

protected:
  CachedNdft<Real, RDmn, WTpExtDmn, WTpExtPosDmn, linalg::CPU, non_density_density_> ndft_obj_;

  int sign_;
private:
  // work spaces for computeGMultiband.
  Matrix G0_M_, G_a_, G_b_;
};

template <class Parameters, DistType DT>
TpAccumulator<Parameters, DT, linalg::CPU>::TpAccumulator(
    const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
    const Parameters& pars, const int thread_id)
    : Base(G0, pars, thread_id),
      G0_M_(n_bands_),
      G_a_(n_bands_),
      G_b_(n_bands_) {
  if constexpr (DT == DistType::BLOCKED) {
    std::cerr << "Blocked distribution is not supported in the CPU accumulator. "
              << "Reverting to no distribution.\n";
  }

}

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::CPU>::resetAccumulation(unsigned int /*dca_loop*/) {
  for (auto& G4_channel : G4_)
    G4_channel = 0.;

  Base::initializeG0();
}

template <class Parameters, DistType DT>
template <class Configuration, typename RealIn>
double TpAccumulator<Parameters, DT, linalg::CPU>::accumulate(
    const std::array<linalg::Matrix<RealIn, linalg::CPU>, 2>& M_pair,
    const std::array<Configuration, 2>& configs, const int sign) {
  Profiler profiler("accumulate", "tp-accumulation", __LINE__, thread_id_);
  double gflops(0.);
  if (!(configs[0].size() + configs[1].size()))  // empty config
    return gflops;

  sign_ = sign;
  gflops += computeM(M_pair, configs);
  gflops += computeG();

  for (int channel_id = 0; channel_id < G4_.size(); ++channel_id)
    gflops += updateG4(channel_id);

  return gflops;
}

template <class Parameters, DistType DT>
template <class Configuration, typename RealIn>
float TpAccumulator<Parameters, DT, linalg::CPU>::computeM(
    const std::array<linalg::Matrix<RealIn, linalg::CPU>, 2>& M_pair,
    const std::array<Configuration, 2>& configs) {
  float flops = 0.;

  func::function<Complex, func::dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, WTpExtPosDmn, WTpExtDmn>> M_r_r_w_w;

  for (int spin = 0; spin < SDmn::dmn_size(); ++spin) {
    Profiler prf_a("Frequency FT", "tp-accumulation", __LINE__, thread_id_);
    if (not configs[spin].size())
      continue;
    flops += ndft_obj_.execute(configs[spin], M_pair[spin], M_r_r_w_w, spin);
  }

  Profiler prf_b("Space FT", "tp-accumulation", __LINE__, thread_id_);
  // TODO: add the gflops here.
  math::transform::SpaceTransform2D<RDmn, KDmn, Real>::execute(M_r_r_w_w, G_);

  return flops;
}

template <class Parameters, DistType DT>
double TpAccumulator<Parameters, DT, linalg::CPU>::computeG() {
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

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::CPU>::computeGSingleband(const int s, const int k1,
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

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::CPU>::computeGMultiband(const int s, const int k1,
                                                                   const int k2, const int w1,
                                                                   const int w2) {
  assert(w1 < WTpExtPosDmn::dmn_size());
  assert(w2 < WTpExtDmn::dmn_size());

  const linalg::MatrixView<Complex, linalg::CPU> G0_w1(&G0_(0, 0, s, k1, w1 + Base::n_pos_frqs_),
                                                       Base::n_bands_, Base::n_bands_);
  const linalg::MatrixView<Complex, linalg::CPU> G0_w2(&G0_(0, 0, s, k2, w2), Base::n_bands_, Base::n_bands_);
  linalg::MatrixView<Complex, linalg::CPU> M_matrix(&G_(0, 0, s, k1, k2, w1, w2), Base::n_bands_);

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

template <class Parameters, DistType DT>
auto TpAccumulator<Parameters, DT, linalg::CPU>::getGSingleband(const int s, const int k1,
                                                                const int k2, const int w1,
                                                                const int w2) -> Complex const {
  const int w2_ext = w2 + extension_index_offset_;
  const int w1_ext = w1 + extension_index_offset_;
  auto minus_w1 = [=](const int w) { return n_pos_frqs_ - 1 - w; };
  auto minus_w2 = [=](const int w) { return 2 * n_pos_frqs_ - 1 - w; };
  auto plus_w1 = [=](const int w) { return w - n_pos_frqs_; };
  auto minus_k = [=](const int k) {
    const static int k0 = KDmn::parameter_type::origin_index();
    return KDmn::parameter_type::subtract(k, k0);
  };

  if (w1_ext >= n_pos_frqs_) {
    // \todo This check can be probably be dropped it worked around another bug
    if (w2_ext < 0)
      return G_(0, 0, s, k1, k2, plus_w1(w1_ext), plus_w1(w2_ext));
    else
      return G_(0, 0, s, k1, k2, plus_w1(w1_ext), w2_ext);
  }
  else
    return std::conj(G_(0, 0, s, minus_k(k1), minus_k(k2), minus_w1(w1_ext), minus_w2(w2_ext)));
}

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::CPU>::getGMultiband(int s, int k1, int k2, int w1, int w2,
                                                               Matrix& G, const Complex beta) const {
  const int w2_ext = w2 + extension_index_offset_;
  const int w1_ext = w1 + extension_index_offset_;
  auto minus_w1 = [=](const int w) { return n_pos_frqs_ - 1 - w; };
  auto minus_w2 = [=](const int w) { return 2 * n_pos_frqs_ - 1 - w; };
  auto plus_w1 = [=](const int w) { return w - n_pos_frqs_; };

  auto minus_k = [=](const int k) {
    const static int k0 = KDmn::parameter_type::origin_index();
    return KDmn::parameter_type::subtract(k, k0);
  };

  if (w1_ext >= n_pos_frqs_) {
    if (w2_ext < 0) {
      const Complex* const G_ptr = &G_(0, 0, s, k1, k2, plus_w1(w1_ext), plus_w1(w2_ext));
      for (int b2 = 0; b2 < n_bands_; ++b2)
        for (int b1 = 0; b1 < n_bands_; ++b1)
          G(b1, b2) = beta * G(b1, b2) + G_ptr[b1 + b2 * n_bands_];
    }
    else {
      const Complex* const G_ptr = &G_(0, 0, s, k1, k2, plus_w1(w1_ext), w2_ext);
      for (int b2 = 0; b2 < n_bands_; ++b2)
        for (int b1 = 0; b1 < n_bands_; ++b1)
          G(b1, b2) = beta * G(b1, b2) + G_ptr[b1 + b2 * n_bands_];
    }
  }
  else {
    const Complex* const G_ptr =
        &G_(0, 0, s, minus_k(k1), minus_k(k2), minus_w1(w1_ext), minus_w2(w2_ext));
    for (int b2 = 0; b2 < n_bands_; ++b2)
      for (int b1 = 0; b1 < n_bands_; ++b1)
        G(b1, b2) = beta * G(b1, b2) + std::conj(G_ptr[b1 + b2 * n_bands_]);
  }
}

template <class Parameters, DistType DT>
double TpAccumulator<Parameters, DT, linalg::CPU>::updateG4(const int channel_id) {
  // G4 is stored with the following band convention:
  // b1 ------------------------ b3
  //        |           |
  //        |           |
  //        |           |
  // b2 ------------------------ b4
  Profiler profiler("updateG4", "tp-accumulation", __LINE__, thread_id_);

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

  auto& G4 = G4_[channel_id];
  auto channel = channels_[channel_id];

  switch (channel) {
    case FourPointType::PARTICLE_HOLE_TRANSVERSE:
      // G4(k1, k2, k_ex) = 1/2 sum_s <c^+(k1+k_ex, s) c(k1, -s) c^+(k2, -s) c(k2+k_ex, s)>
      //                  = -1/2 sum_s G(k2+k_ex, k1+k_ex, s) G(k1, k2, -s)
      for (int w_ex_idx = 0; w_ex_idx < exchange_frq.size(); ++w_ex_idx) {
        const int w_ex = exchange_frq[w_ex_idx];
        for (int k_ex_idx = 0; k_ex_idx < exchange_mom.size(); ++k_ex_idx) {
          const int k_ex = exchange_mom[k_ex_idx];
          for (int w2 = 0; w2 < WTpDmn::dmn_size(); ++w2)
            for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
              for (int w1 = 0; w1 < WTpDmn::dmn_size(); ++w1)
                for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
                  Complex* const G4_ptr = &G4(0, 0, 0, 0, k1, w1, k2, w2, k_ex_idx, w_ex_idx);
                  for (int s = 0; s < 2; ++s)
                    updateG4Atomic(G4_ptr, s, k1, k2, w1, w2, not s, momentum_sum(k2, k_ex),
                                   momentum_sum(k1, k_ex), w_plus_w_ex(w2, w_ex),
                                   w_plus_w_ex(w1, w_ex), -sign_over_2, true);
                }
        }
      }
      flops += n_loops * 2 * flops_update_atomic;
      break;

    case FourPointType::PARTICLE_HOLE_MAGNETIC:
      // G4(k1, k2, k_ex) = 1/2 sum_{s1, s2} (s1 * s2)
      //                      <c^+(k1+k_ex, s1) c(k1, s1) c^+(k2, s2) c(k2+k_ex, s2)>
      //                  = 1/2 sum_{s1, s2} (s1 * s2)
      //                      [G(k1, k1+k_ex, s1) G(k2+k_ex, k2, s2)
      //                       - (s1 == s2) G(k2+k_ex, k1+k_ex, s1) G(k1, k2, s1)]
      for (int w_ex_idx = 0; w_ex_idx < exchange_frq.size(); ++w_ex_idx) {
        const int w_ex = exchange_frq[w_ex_idx];
        for (int k_ex_idx = 0; k_ex_idx < exchange_mom.size(); ++k_ex_idx) {
          const int k_ex = exchange_mom[k_ex_idx];
          for (int w2 = 0; w2 < WTpDmn::dmn_size(); ++w2)
            for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
              for (int w1 = 0; w1 < WTpDmn::dmn_size(); ++w1)
                for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
                  Complex* const G4_ptr = &G4(0, 0, 0, 0, k1, w1, k2, w2, k_ex_idx, w_ex_idx);
                  updateG4SpinDifference(G4_ptr, -1, k1, momentum_sum(k1, k_ex), w1,
                                         w_plus_w_ex(w1, w_ex), momentum_sum(k2, k_ex), k2,
                                         w_plus_w_ex(w2, w_ex), w2, sign_over_2, false);
                  for (int s = 0; s < 2; ++s)
                    updateG4Atomic(G4_ptr, s, k1, k2, w1, w2, s, momentum_sum(k2, k_ex),
                                   momentum_sum(k1, k_ex), w_plus_w_ex(w2, w_ex),
                                   // w_plus_w_ex(w1, w_ex), -sign_over_2, true);
                                   w_plus_w_ex(w1, w_ex), -sign_over_2, false);
                }
        }
      }
      flops += n_loops * (flops_update_spin_diff + 2 * flops_update_atomic);
      break;

    case FourPointType::PARTICLE_HOLE_CHARGE:
      // G4(k1, k2, k_ex) = 1/2 sum_{s1, s2}
      //                    <c^+(k1+k_ex, s1) c(k1, s1) c^+(k2, s2) c(k2+k_ex, s2)>
      //                  = 1/2 sum_{s1, s2}
      //                      [G(k1, k1+k_ex, s1) G(k2+k_ex, k2, s2)
      //                       - (s1 == s2) G(k2+k_ex, k1+k_ex, s1) G(k1, k2, s1)]
      for (int w_ex_idx = 0; w_ex_idx < exchange_frq.size(); ++w_ex_idx) {
        const int w_ex = exchange_frq[w_ex_idx];
        for (int k_ex_idx = 0; k_ex_idx < exchange_mom.size(); ++k_ex_idx) {
          const int k_ex = exchange_mom[k_ex_idx];
          for (int w2 = 0; w2 < WTpDmn::dmn_size(); ++w2)
            for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
              for (int w1 = 0; w1 < WTpDmn::dmn_size(); ++w1)
                for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
                  Complex* const G4_ptr = &G4(0, 0, 0, 0, k1, w1, k2, w2, k_ex_idx, w_ex_idx);
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

    case FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP:
      // G4(k1, k2, k_ex) = 1/2 sum_s <c^+(k1+k_ex, s) c(k1, s) c^+(k2, s) c(k2+k_ex, s)>
      //                  = 1/2 sum_s [G(k1, k1+k_ex, s) G(k2+k_ex, k2, s)
      //                               - G(k2+k_ex, k1+k_ex, s) G(k1, k2, s)]
      for (int w_ex_idx = 0; w_ex_idx < exchange_frq.size(); ++w_ex_idx) {
        const int w_ex = exchange_frq[w_ex_idx];
        for (int k_ex_idx = 0; k_ex_idx < exchange_mom.size(); ++k_ex_idx) {
          const int k_ex = exchange_mom[k_ex_idx];
          for (int w2 = 0; w2 < WTpDmn::dmn_size(); ++w2)
            for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
              for (int w1 = 0; w1 < WTpDmn::dmn_size(); ++w1)
                for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
                  Complex* const G4_ptr = &G4(0, 0, 0, 0, k1, w1, k2, w2, k_ex_idx, w_ex_idx);

                  for (int s = 0; s < 2; ++s)
                    updateG4Atomic(G4_ptr, s, k1, momentum_sum(k1, k_ex), w1, w_plus_w_ex(w1, w_ex),
                                   s, momentum_sum(k2, k_ex), k2, w_plus_w_ex(w2, w_ex), w2,
                                   sign_over_2, false);

                  for (int s = 0; s < 2; ++s)
                    updateG4Atomic(G4_ptr, s, k1, k2, w1, w2, s, momentum_sum(k2, k_ex),
                                   momentum_sum(k1, k_ex), w_plus_w_ex(w2, w_ex),
                                   w_plus_w_ex(w1, w_ex), -sign_over_2, true);
                }
        }
      }
      flops += n_loops * 4 * flops_update_atomic;
      break;

    case FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN:
      // G4(k1, k2, k_ex) = 1/2 sum_s <c^+(k1+k_ex, s) c(k1, s) c^+(k2, -s) c(k2+k_ex, -s)>
      //                  = 1/2 sum_s G(k1, k1+k_ex, s) G(k2+k_ex, k2, -s)
      for (int w_ex_idx = 0; w_ex_idx < exchange_frq.size(); ++w_ex_idx) {
        const int w_ex = exchange_frq[w_ex_idx];
        for (int k_ex_idx = 0; k_ex_idx < exchange_mom.size(); ++k_ex_idx) {
          const int k_ex = exchange_mom[k_ex_idx];
          for (int w2 = 0; w2 < WTpDmn::dmn_size(); ++w2)
            for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
              for (int w1 = 0; w1 < WTpDmn::dmn_size(); ++w1)
                for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
                  Complex* const G4_ptr = &G4(0, 0, 0, 0, k1, w1, k2, w2, k_ex_idx, w_ex_idx);
                  for (int s = 0; s < 2; ++s)
                    updateG4Atomic(G4_ptr, s, k1, momentum_sum(k1, k_ex), w1, w_plus_w_ex(w1, w_ex),
                                   !s, momentum_sum(k2, k_ex), k2, w_plus_w_ex(w2, w_ex), w2,
                                   sign_over_2, false);
                }
        }
      }
      flops += n_loops * 4 * flops_update_atomic;
      break;

    case FourPointType::PARTICLE_PARTICLE_UP_DOWN:
      // G4(k1, k2, k_ex) = 1/2 sum_s <c^+(k_ex-k1, s) c^+(k1, -s) c(k2, -s) c(k_ex-k2, s)>
      //                  = 1/2 sum_s G(k_ex-k2, k_ex-k1, s) G(k2, k1, -s)
      for (int w_ex_idx = 0; w_ex_idx < exchange_frq.size(); ++w_ex_idx) {
        const int w_ex = exchange_frq[w_ex_idx];
        for (int k_ex_idx = 0; k_ex_idx < exchange_mom.size(); ++k_ex_idx) {
          const int k_ex = exchange_mom[k_ex_idx];
          for (int w2 = 0; w2 < WTpDmn::dmn_size(); ++w2)
            for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
              for (int w1 = 0; w1 < WTpDmn::dmn_size(); ++w1)
                for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
                  Complex* const G4_ptr = &G4(0, 0, 0, 0, k1, w1, k2, w2, k_ex_idx, w_ex_idx);
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
      throw std::logic_error("Specified four point type not implemented.");
  }

  return 1e-9 * flops;
}

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::CPU>::updateG4Atomic(
    Complex* G4_ptr, const int s_a, const int k1_a, const int k2_a, const int w1_a, const int w2_a,
    const int s_b, const int k1_b, const int k2_b, const int w1_b, const int w2_b, const Real alpha,
    const bool cross_legs) {
  // This function performs the following update for each band:
  //
  // G4(k1, k2, w1, w2) += alpha * G(s_a, k1_a, k2_a, w1_a, w2_a) * G(s_b, k1_b, k2_b, w1_b, w2_b)
  //
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
              // *G4_ptr += alpha * G_a_(b1, b3) * G_b_(b2, b4);
              *G4_ptr += alpha * G_a_(b2, b4) * G_b_(b3, b1);
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

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::CPU>::updateG4SpinDifference(
    Complex* G4_ptr, const int sign, const int k1_a, const int k2_a, const int w1_a, const int w2_a,
    const int k1_b, const int k2_b, const int w1_b, const int w2_b, const Real alpha,
    const bool cross_legs) {
  // This function performs the following update for each band:
  //
  // G4(k1, k2, w1, w2) += alpha * [G(up, k1_a, k2_a, w1_a, w2_a)
  //                                + sign * G(down,k1_a, k2_a, w1_a, w2_a)]
  //                             * [G(up, k1_b, k2_b, w1_b, w2_b)
  //                               + sign * G(down, k1_b, k2_b, w1_b, w2_b)]
  if (n_bands_ == 1) {
    *G4_ptr += alpha *
               (getGSingleband(0, k1_a, k2_a, w1_a, w2_a) +
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
              // *G4_ptr += alpha * G_a_(b1, b3) * G_b_(b2, b4);
              *G4_ptr += alpha * G_a_(b2, b1) * G_b_(b3, b4);
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

template <class Parameters, DistType DT>
const std::vector<typename TpAccumulator<Parameters, DT, linalg::CPU>::TpGreensFunction>& TpAccumulator<
    Parameters, DT, linalg::CPU>::get_G4() const {
  if (G4_.empty())
    throw std::logic_error("There is no G4 stored in this class.");
  return G4_;
}

template <class Parameters, DistType DT>
std::vector<typename TpAccumulator<Parameters, DT, linalg::CPU>::TpGreensFunction>& TpAccumulator<
    Parameters, DT, linalg::CPU>::get_nonconst_G4() {
  if (G4_.empty())
    throw std::logic_error("There is no G4 stored in this class.");
  return G4_;
}

template <class Parameters, DistType DT>
void TpAccumulator<Parameters, DT, linalg::CPU>::sumTo(TpAccumulator& other_one) {
  if (other_one.G4_.size() != G4_.size())
    throw std::logic_error("Objects accumulate different number of channels.");

  for (std::size_t channel = 0; channel < G4_.size(); ++channel)
    other_one.G4_[channel] += G4_[channel];

  G4_.clear();
}

}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_HPP
