// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This class provides a base for the device and G4 distribution specializations
// that compute the two particle Green's function from the walker's M matrix.

#ifndef DCA_TP_ACCUMULATOR_BASE_HPP
#define DCA_TP_ACCUMULATOR_BASE_HPP

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

using dca::util::RealAlias;
  
template <class Parameters, DistType DT = DistType::NONE, linalg::DeviceType device = linalg::CPU>
class TpAccumulator {};

template <class Parameters, DistType DT>
class TpAccumulatorBase {
public:
  using Real = typename Parameters::TP_measurement_scalar_type;
  using Scalar = typename dca::util::ScalarSelect<Real,Parameters::complex_g0>::type;
  
  using RDmn = typename Parameters::RClusterDmn;
  using KDmn = typename Parameters::KClusterDmn;
  using KExchangeDmn = func::dmn_0<domains::MomentumExchangeDomain>;
  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;
  using WDmn = func::dmn_0<domains::frequency_domain>;

  using WTpDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using WTpPosDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT_POSITIVE>>;
  using WTpExtDmn = func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED>>;
  using WTpExtPosDmn = func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED_POSITIVE>>;
  using WExchangeDmn = func::dmn_0<domains::FrequencyExchangeDomain>;

  using TpGreensFunction = typename DcaData<Parameters, DT>::TpGreensFunction;
  
protected:
  using Profiler = typename Parameters::profiler_type;

  using Complex = std::complex<RealAlias<Scalar>>;

  using SpGreenFunction =
      func::function<Complex, func::dmn_variadic<BDmn, BDmn, SDmn, KDmn, KDmn, WTpExtPosDmn, WTpExtDmn>>;

  using TpDomain =
      func::dmn_variadic<BDmn, BDmn, BDmn, BDmn, KDmn, KDmn, KExchangeDmn, WTpDmn, WTpDmn, WExchangeDmn>;

public:
  // Constructor:
  // In: G0: non interacting greens function.
  // In: pars: parameters object.
  // In: thread_id: thread id, only used by the profiler.
  TpAccumulatorBase(
      const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
      const Parameters& pars, int thread_id = 0);

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

  // Sums the accumulated Green's function to the accumulated Green's function of other_acc.
  void sumTo(TpAccumulatorBase<Parameters, DT>& other_acc);

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

  auto get_n_pos_frqs() {
    return n_pos_frqs_;
  }
protected:
  void initializeG0();

  double computeG();

  void computeGMultiband(int s, int k1, int k2, int w1, int w2);

  void computeGSingleband(int s, int k1, int k2, int w1, int w2);

//  void getGMultiband(int s, int k1, int k2, int w1, int w2, Matrix& G, Complex beta = 0) const;

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
  const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>* const G0_ptr_ =
      nullptr;

  const int thread_id_;
  bool multiple_accumulators_;

  const Real beta_ = -1;
  constexpr static int n_bands_ = Parameters::model_type::BANDS;

  constexpr static bool  spin_symmetric_ = Parameters::model_type::spin_symmetric;
  
  constexpr static bool non_density_density_ =
    models::HasInitializeNonDensityInteractionMethod<Parameters>::value;
  CachedNdft<Real, RDmn, WTpExtDmn, WTpExtPosDmn, linalg::CPU, non_density_density_> ndft_obj_;

  SpGreenFunction G_;

  std::vector<TpGreensFunction> G4_;
  std::vector<FourPointType> channels_;

  func::function<Complex, func::dmn_variadic<BDmn, BDmn, SDmn, KDmn, WTpExtDmn>> G0_;

  math::Phase<Scalar> phase_;

  const int extension_index_offset_ = -1;
  const int n_pos_frqs_ = -1;

};

template <class Parameters, DistType DT>
TpAccumulatorBase<Parameters, DT>::TpAccumulatorBase(
    const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
    const Parameters& pars, const int thread_id)
    : G0_ptr_(&G0),
      thread_id_(thread_id),
      multiple_accumulators_(pars.get_accumulators() > 1),
      beta_(pars.get_beta()),
      channels_(pars.get_four_point_channels()),
      extension_index_offset_((WTpExtDmn::dmn_size() - WTpDmn::dmn_size()) / 2),
      n_pos_frqs_(WTpExtPosDmn::dmn_size()) {

  if (WDmn::dmn_size() < WTpExtDmn::dmn_size())
    throw(std::logic_error("The number of single particle frequencies is too small."));

  initializeG0();

  // Reserve storage in advance such that we don't have to copy elements when we fill the vector.
  // We want to avoid copies because function's copy ctor does not copy the name (and because copies
  // are expensive).
  for (auto channel : channels_) {
    G4_.emplace_back("G4_" + toString(channel), pars.get_concurrency());
  }

}

template <class Parameters, DistType DT>
void TpAccumulatorBase<Parameters, DT>::initializeG0() {
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

}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_TP_ACCUMULATOR_BASH_HPP
