// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements a 2d NDFT while tracking the band indices of each measurement.

#ifndef DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_HPP
#define DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_HPP

#include <array>
#include <cassert>
#include <complex>
#include <type_traits>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/linalg/util/allocators.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/triple.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/util/ignore.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, linalg::DeviceType device,
          bool non_density_density = 0>
class CachedNdft;

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
class CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density> {
private:
  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using ClusterDmn = typename RDmn::parameter_type;
  using BRDmn = func::dmn_variadic<BDmn, RDmn>;
  using Matrix = linalg::Matrix<ScalarType, dca::linalg::CPU>;

public:
  CachedNdft();

  // For each pair of orbitals, performs the non-uniform 2D Fourier Transform from time to frequency
  // defined as M(w1, w2) = \sum_{t1, t2} exp(i (w1 t1 - w2 t2)) M(t1, t2).
  // In: configuration: stores the time and orbitals of each entry of M.
  // In: M: input matrix provided by the walker.
  // Out: M_r_r_w_w: stores the result of the computation.
  // In: spin: spin sector of M. If required it is used to sort the output.
  template <class Configuration, typename ScalarInp, class OutDmn>
  double execute(const Configuration& configuration, const linalg::Matrix<ScalarInp, linalg::CPU>& M,
                 func::function<std::complex<ScalarType>, OutDmn>& M_r_r_w_w, int spin = 0);

protected:
  template <class Configuration>
  void sortConfiguration(const Configuration& configuration);

private:
  template <class Configuration>
  void computeT(const Configuration& configuration);

  template <class Configuration, typename ScalarInp>
  void computeMMatrix(const Configuration& configuration,
                      const linalg::Matrix<ScalarInp, linalg::CPU>& M, int orb_i, int orb_j);

  void computeTSubmatrices(int orb_i, int orb_j);

  double executeTrimmedFT();

  void inline copyPartialResult(
      int orb1, int orb2, int /*spin*/,
      func::function<std::complex<ScalarType>,
                     func::dmn_variadic<BDmn, BDmn, RDmn, RDmn, WPosDmn, WDmn>>& f_out) const;

  void inline copyPartialResult(
      int orb1, int orb2, int spin,
      func::function<std::complex<ScalarType>,
                     func::dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, WPosDmn, WDmn>>& f_out) const;

  void inline setToZero(
      int orb1, int orb2, int /*spin*/,
      func::function<std::complex<ScalarType>,
                     func::dmn_variadic<BDmn, BDmn, RDmn, RDmn, WPosDmn, WDmn>>& f_out) const;

  void inline setToZero(
      int orb1, int orb2, int spin,
      func::function<std::complex<ScalarType>,
                     func::dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, WPosDmn, WDmn>>& f_out) const;

  static void orbitalToBR(int orbital, int& b, int& r);

protected:
  using Triple = details::Triple<ScalarType>;
  template <class T>
  using HostVector = linalg::util::HostVector<T>;

  HostVector<ScalarType> w_;

  std::array<HostVector<Triple>, 2> indexed_config_;

  std::array<std::vector<int>, 2> start_index_;
  std::array<std::vector<int>, 2> end_index_;

  const HostVector<Triple> &config_left_, &config_right_;
  std::vector<int> &start_index_left_, &start_index_right_;
  std::vector<int> &end_index_left_, &end_index_right_;

  const int n_orbitals_;

private:
  using MatrixPair = std::array<Matrix, 2>;
  Matrix M_ij_;
  MatrixPair T_l_times_M_ij_times_T_r_;
  MatrixPair T_;
  MatrixPair T_l_;
  MatrixPair T_r_;
  MatrixPair T_l_times_M_ij_;
  std::array<Matrix, 5> work_;
};

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::CachedNdft()
    : w_(),
      config_left_(indexed_config_[0]),
      config_right_(non_density_density ? indexed_config_[1] : indexed_config_[0]),
      start_index_left_(start_index_[0]),
      start_index_right_(non_density_density ? start_index_[1] : start_index_[0]),
      end_index_left_(end_index_[0]),
      end_index_right_(non_density_density ? end_index_[1] : end_index_[0]),
      n_orbitals_(BDmn::dmn_size() * RDmn::dmn_size()) {
  for (const auto elem : WDmn::parameter_type::get_elements())
    w_.push_back(static_cast<ScalarType>(elem));

  start_index_left_.resize(n_orbitals_);
  start_index_right_.resize(n_orbitals_);
  end_index_right_.resize(n_orbitals_);
  end_index_right_.resize(n_orbitals_);
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
template <class Configuration, typename ScalarInp, class OutDmn>
double CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::execute(
    const Configuration& configuration, const linalg::Matrix<ScalarInp, linalg::CPU>& M,
    func::function<std::complex<ScalarType>, OutDmn>& M_r_r_w_w, const int spin) {
  assert(M_r_r_w_w[M_r_r_w_w.signature() - 1] == WDmn::dmn_size());
  assert(M_r_r_w_w[M_r_r_w_w.signature() - 2] == WPosDmn::dmn_size());

  double flops = 0.;
  //  const int r0_index = ClusterDmn::origin_index();

  sortConfiguration(configuration);

  computeT(configuration);

  for (int orb_j = 0; orb_j < n_orbitals_; orb_j++) {
    //  INTERNAL: Is there any case where r_j - r_0 != r_j?
    //  int min_r_j = ClusterDmn::subtract(r_j, r0_index);
    const int n_j = end_index_right_[orb_j] - start_index_right_[orb_j];

    for (int orb_i = 0; orb_i < n_orbitals_; orb_i++) {
      const int n_i = end_index_left_[orb_i] - start_index_left_[orb_i];

      if (n_i > 0 && n_j > 0) {
        computeMMatrix(configuration, M, orb_i, orb_j);

        computeTSubmatrices(orb_i, orb_j);

        flops += executeTrimmedFT();

        copyPartialResult(orb_i, orb_j, spin, M_r_r_w_w);
      }
      else
        setToZero(orb_i, orb_j, spin, M_r_r_w_w);
    }
  }

  return (flops * (1.e-9));
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
template <class Configuration>
void CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::sortConfiguration(
    const Configuration& configuration) {
  const int n_b = BDmn::dmn_size();
  const int n_v = configuration.size();
  constexpr int n_sides = non_density_density ? 2 : 1;

  auto orbital = [&](const int i, const int side) {
    const auto& vertex = configuration[i];
    if (side)  // Switch left and right band as M is an inverse matrix.
      return vertex.get_left_band() + n_b * vertex.get_left_site();
    else
      return vertex.get_right_band() + n_b * vertex.get_right_site();
  };

  for (int side = 0; side < n_sides; ++side) {
    auto& config_side = indexed_config_[side];
    config_side.resize(n_v);

    for (int l = 0; l < n_v; l++) {
      config_side[l].orbital = orbital(l, side);
      config_side[l].tau = configuration[l].get_tau();
      config_side[l].idx = l;
    }

    sort(config_side.begin(), config_side.end());

    for (int orb = 0; orb < n_orbitals_; orb++) {
      details::Triple<ScalarType> trp{orb, 0, 0};

      start_index_[side][orb] =
          lower_bound(config_side.begin(), config_side.end(), trp) - config_side.begin();

      trp.idx = n_v;

      end_index_[side][orb] =
          upper_bound(config_side.begin(), config_side.end(), trp) - config_side.begin();
    }
  }
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
template <class Configuration>
void CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::computeT(
    const Configuration& configuration) {
  int n_v = configuration.size();
  int n_w = w_.size();

  T_[0].resizeNoCopy(std::pair<int, int>(n_w, n_v));
  T_[1].resizeNoCopy(std::pair<int, int>(n_w, n_v));

  for (int j = 0; j < n_v; j++) {
    for (int i = 0; i < n_w; i++) {
      const ScalarType x = configuration[j].get_tau() * w_[i];

      T_[0](i, j) = cos(x);
      T_[1](i, j) = sin(x);
    }
  }
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
template <class Configuration, typename ScalarInp>
void CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::computeMMatrix(
    const Configuration& configuration, const linalg::Matrix<ScalarInp, linalg::CPU>& M,
    const int orb_i, const int orb_j) {
  // In release mode 'configuration' is an unused parameter.
  dca::util::ignoreUnused(configuration);

  M_ij_.resizeNoCopy(std::pair<int, int>(end_index_left_[orb_i] - start_index_left_[orb_i],
                                         end_index_right_[orb_j] - start_index_right_[orb_j]));

  for (int l_j = start_index_right_[orb_j]; l_j < end_index_right_[orb_j]; ++l_j) {
    const int out_j = l_j - start_index_right_[orb_j];
    for (int l_i = start_index_left_[orb_i]; l_i < end_index_left_[orb_i]; ++l_i) {
      const int out_i = l_i - start_index_left_[orb_i];
      M_ij_(out_i, out_j) = M(config_left_[l_i].idx, config_right_[l_j].idx);
    }
  }
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
void CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::computeTSubmatrices(
    const int orb_i, const int orb_j) {
  const int n_w = WDmn::dmn_size();
  const int n_w_pos = WPosDmn::dmn_size();

  for (int re_im = 0; re_im < 2; ++re_im) {
    // T_l_ matrix
    T_l_[re_im].resizeNoCopy(
        std::pair<int, int>(n_w_pos, end_index_left_[orb_i] - start_index_left_[orb_i]));

    for (int l_i = start_index_left_[orb_i]; l_i < end_index_left_[orb_i]; ++l_i) {
      const int i = l_i - start_index_left_[orb_i];
      memcpy(&T_l_[re_im](0, i), &T_[re_im](n_w_pos, config_left_[l_i].idx),
             sizeof(ScalarType) * n_w_pos);
    }

    // T_r_ matrix
    T_r_[re_im].resizeNoCopy(
        std::pair<int, int>(n_w, end_index_right_[orb_j] - start_index_right_[orb_j]));

    for (int l_j = start_index_right_[orb_j]; l_j < end_index_right_[orb_j]; ++l_j) {
      const int j = l_j - start_index_right_[orb_j];
      memcpy(&T_r_[re_im](0, j), &T_[re_im](0, config_right_[l_j].idx), sizeof(ScalarType) * n_w);
    }
  }
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
double CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::executeTrimmedFT() {
  double flops = 0.;

  assert(WPosDmn::dmn_size() == WDmn::dmn_size() / 2);

  assert(T_l_[0].size().first == WPosDmn::dmn_size());
  assert(T_l_[0].size().second == M_ij_.size().first);

  assert(T_r_[0].size().first == WDmn::dmn_size());
  assert(T_r_[0].size().second == M_ij_.size().second);

  for (int re_im = 0; re_im < 2; ++re_im) {
    T_l_times_M_ij_[re_im].resizeNoCopy(std::make_pair(WPosDmn::dmn_size(), M_ij_.size().second));
    T_l_times_M_ij_times_T_r_[re_im].resizeNoCopy(
        std::make_pair(WPosDmn::dmn_size(), WDmn::dmn_size()));
  }

  dca::linalg::matrixop::multiply(T_l_, M_ij_, T_l_times_M_ij_);
  flops += 4 * T_l_[0].size().first * T_l_[0].size().second * M_ij_.size().second;

  dca::linalg::matrixop::multiply('N', 'C', T_l_times_M_ij_, T_r_, T_l_times_M_ij_times_T_r_, work_);
  flops += 4. * T_l_times_M_ij_[0].size().first * T_l_times_M_ij_[0].size().second *
           T_l_times_M_ij_times_T_r_[0].size().second;

  return 1e-9 * flops;
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
void CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::copyPartialResult(
    const int orb1, const int orb2, int /*spin*/,
    func::function<std::complex<ScalarType>,
                   func::dmn_variadic<BDmn, BDmn, RDmn, RDmn, WPosDmn, WDmn>>& f_out) const {
  const int n_w1 = WPosDmn::dmn_size();
  const int n_w2 = WDmn::dmn_size();
  int b1, b2, r1, r2;
  orbitalToBR(orb1, b1, r1);
  orbitalToBR(orb2, b2, r2);
  for (int w2 = 0; w2 < n_w2; ++w2)
    for (int w1 = 0; w1 < n_w1; ++w1) {
      f_out(b1, b2, r1, r2, w1, w2).real(T_l_times_M_ij_times_T_r_[0](w1, w2));
      f_out(b1, b2, r1, r2, w1, w2).imag(T_l_times_M_ij_times_T_r_[1](w1, w2));
    }
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
void CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::copyPartialResult(
    const int orb1, const int orb2, const int spin,
    func::function<std::complex<ScalarType>,
                   func::dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, WPosDmn, WDmn>>& f_out) const {
  const int n_w1 = WPosDmn::dmn_size();
  const int n_w2 = WDmn::dmn_size();
  int b1, b2, r1, r2;
  orbitalToBR(orb1, b1, r1);
  orbitalToBR(orb2, b2, r2);
  for (int w2 = 0; w2 < n_w2; ++w2)
    for (int w1 = 0; w1 < n_w1; ++w1) {
      f_out(r1, r2, b1, b2, spin, w1, w2).real(T_l_times_M_ij_times_T_r_[0](w1, w2));
      f_out(r1, r2, b1, b2, spin, w1, w2).imag(T_l_times_M_ij_times_T_r_[1](w1, w2));
    }
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
void CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::setToZero(
    const int orb1, const int orb2, int /*spin*/,
    func::function<std::complex<ScalarType>,
                   func::dmn_variadic<BDmn, BDmn, RDmn, RDmn, WPosDmn, WDmn>>& f_out) const {
  const int n_w1 = WPosDmn::dmn_size();
  const int n_w2 = WDmn::dmn_size();
  int b1, b2, r1, r2;
  orbitalToBR(orb1, b1, r1);
  orbitalToBR(orb2, b2, r2);
  for (int w2 = 0; w2 < n_w2; ++w2)
    for (int w1 = 0; w1 < n_w1; ++w1)
      f_out(b1, b2, r1, r2, w1, w2) = 0;
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
void CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::setToZero(
    const int orb1, const int orb2, const int spin,
    func::function<std::complex<ScalarType>,
                   func::dmn_variadic<RDmn, RDmn, BDmn, BDmn, SDmn, WPosDmn, WDmn>>& f_out) const {
  const int n_w1 = WPosDmn::dmn_size();
  const int n_w2 = WDmn::dmn_size();
  int b1, b2, r1, r2;
  orbitalToBR(orb1, b1, r1);
  orbitalToBR(orb2, b2, r2);
  for (int w2 = 0; w2 < n_w2; ++w2)
    for (int w1 = 0; w1 < n_w1; ++w1)
      f_out(r1, r2, b1, b2, spin, w1, w2) = 0;
}

template <typename ScalarType, class RDmn, class WDmn, class WPosDmn, bool non_density_density>
void CachedNdft<ScalarType, RDmn, WDmn, WPosDmn, linalg::CPU, non_density_density>::orbitalToBR(
    int orbital, int& b, int& r) {
  const static int n_bands = BDmn::dmn_size();
  r = orbital / n_bands;
  b = orbital % n_bands;
}

}  // accumulator
}  // solver
}  // phys
}  // dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_NDFT_CACHED_NDFT_HPP
