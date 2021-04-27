// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This class provides a base for the device and G4 distribution specializations
// that compute the two particle Green's function from the walker's M matrix.

#ifndef DCA_TP_ACCUMULATOR_GPU_BASE_HPP
#define DCA_TP_ACCUMULATOR_GPU_BASE_HPP

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "dca/config/mc_options.hpp"
#include "dca/distribution/dist_types.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrix_view.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/linalg/util/cuda_stream.hpp"
#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/reshapable_matrix.hpp"
#include "dca/linalg/util/allocators/managed_allocator.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/util/magma_queue.hpp"

#include "dca/math/function_transform/special_transforms/space_transform_2D_gpu.hpp"
#include "dca/phys/dca_data/dca_data.hpp"

#include "dca/phys/domains/cluster/momentum_exchange_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_exchange_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/models/traits.hpp"
#include "dca/phys/four_point_type.hpp"

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/g4_helper.cuh"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {

template <class Parameters, DistType DT>
class TpAccumulatorGpuBase {
public:
  using Real = typename Parameters::TP_measurement_scalar_type;
  using Complex = std::complex<Real>;
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

protected:
  using Profiler = typename Parameters::profiler_type;
  using Matrix = linalg::Matrix<Complex, linalg::GPU>;

  using MatrixDev = linalg::Matrix<Complex, linalg::GPU>;
  using RMatrix =
      linalg::ReshapableMatrix<Complex, linalg::GPU, config::McOptions::TpAllocator<Complex>>;
  using RMatrixValueType = typename RMatrix::ValueType;
  using MatrixHost = linalg::Matrix<Complex, linalg::CPU>;

public:
  TpAccumulatorGpuBase(
      const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
      const Parameters& pars, int n_pos_frqs, int thread_id);

protected:
  void initializeG4Helpers() const;
  void synchronizeStreams();
  void initializeG0();

  template <class Configuration, typename RealIn>
  float computeM(const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M_pair,
                 const std::array<Configuration, 2>& configs);

  void sumTo_(TpAccumulatorGpuBase<Parameters, DT>& other_acc);

  // \todo is this violation of single source of truth necessary.
  const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>* const G0_ptr_ =
      nullptr;

  const int n_pos_frqs_ = -1;

  std::array<linalg::util::MagmaQueue, 2> queues_;
  linalg::util::CudaEvent event_;

  std::vector<std::shared_ptr<RMatrix>> workspaces_;

  // this is how this is defined in the all tp_accumulator_gpu, suspect?
  constexpr static bool non_density_density_ =
      models::has_non_density_interaction<typename Parameters::lattice_type>;
  CachedNdft<Real, RDmn, WTpExtDmn, WTpExtPosDmn, linalg::CPU, non_density_density_> ndft_obj_;

  using NdftType = CachedNdft<Real, RDmn, WTpExtDmn, WTpExtPosDmn, linalg::GPU, non_density_density_>;
  std::array<NdftType, 2> ndft_objs_;
  using DftType = math::transform::SpaceTransform2DGpu<RDmn, KDmn, Real>;
  std::array<DftType, 2> space_trsf_objs_;

  std::array<RMatrix, 2> G_;

  const int nr_accumulators_;

  bool finalized_ = false;
  bool initialized_ = false;

  constexpr static int n_ndft_queues_ = config::McOptions::memory_savings ? 1 : 2;
  constexpr static int n_bands_ = Parameters::model_type::BANDS;

  const int thread_id_;

  using G0DevType = std::array<MatrixDev, 2>;
  static inline G0DevType& get_G0();
};

template <class Parameters, DistType DT>
TpAccumulatorGpuBase<Parameters, DT>::TpAccumulatorGpuBase(
    const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
    const Parameters& pars, const int n_pos_frqs, int thread_id)
    : G0_ptr_(&G0),
      n_pos_frqs_(n_pos_frqs),
      queues_(),
      ndft_objs_{NdftType(queues_[0]), NdftType(queues_[1])},
      space_trsf_objs_{DftType(n_pos_frqs_, queues_[0]), DftType(n_pos_frqs_, queues_[1])},
      nr_accumulators_(pars.get_accumulators()),
      thread_id_(thread_id) {
  initializeG4Helpers();

  // Create shared workspaces.
  for (int i = 0; i < n_ndft_queues_; ++i) {
    workspaces_.emplace_back(std::make_shared<RMatrix>());
    workspaces_[i]->setStream(queues_[i]);
    ndft_objs_[i].setWorkspace(workspaces_[i]);
    space_trsf_objs_[i].setWorkspace(workspaces_[i]);
  }
}

template <class Parameters, DistType DT>
auto TpAccumulatorGpuBase<Parameters, DT>::get_G0() -> G0DevType& {
  static G0DevType G0;
  return G0;
}

template <class Parameters, DistType DT>
void TpAccumulatorGpuBase<Parameters, DT>::initializeG4Helpers() const {
  static std::once_flag flag;
  std::call_once(flag, []() {
    const auto& add_mat = KDmn::parameter_type::get_add_matrix();
    const auto& sub_mat = KDmn::parameter_type::get_subtract_matrix();
    const auto& w_indices = domains::FrequencyExchangeDomain::get_elements();
    const auto& q_indices = domains::MomentumExchangeDomain::get_elements();
    details::G4Helper::set(n_bands_, KDmn::dmn_size(), WTpPosDmn::dmn_size(), q_indices, w_indices,
                           add_mat.ptr(), add_mat.leadingDimension(), sub_mat.ptr(),
                           sub_mat.leadingDimension());
    assert(cudaPeekAtLastError() == cudaSuccess);
  });
}

template <class Parameters, DistType DT>
void TpAccumulatorGpuBase<Parameters, DT>::synchronizeStreams() {
  for (auto& stream : queues_)
    cudaStreamSynchronize(stream);
}

template <class Parameters, DistType DT>
void TpAccumulatorGpuBase<Parameters, DT>::initializeG0() {
  // Note: this method is not thread safe by itself.
  const int sp_index_offset = (WDmn::dmn_size() - WTpExtDmn::dmn_size()) / 2;
  static func::dmn_variadic<BDmn, KDmn, WTpExtDmn> bkw_dmn;
  std::array<MatrixHost, 2> G0_host;

  for (int s = 0; s < 2; ++s) {
    auto& G0 = get_G0();

    G0_host[s].resizeNoCopy(std::make_pair(bkw_dmn.get_size(), BDmn::dmn_size()));
    for (int w = 0; w < WTpExtDmn::dmn_size(); ++w)
      for (int k = 0; k < KDmn::dmn_size(); ++k)
        for (int b2 = 0; b2 < n_bands_; ++b2)
          for (int b1 = 0; b1 < n_bands_; ++b1)
            G0_host[s](bkw_dmn(b1, k, w), b2) = (*G0_ptr_)(b1, s, b2, s, k, w + sp_index_offset);

    G0[s].setAsync(G0_host[s], queues_[s].getStream());
  }
}

template <class Parameters, DistType DT>
template <class Configuration, typename RealIn>
float TpAccumulatorGpuBase<Parameters, DT>::computeM(
    const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M_pair,
    const std::array<Configuration, 2>& configs) {
  auto stream_id = [&](const int s) { return n_ndft_queues_ == 1 ? 0 : s; };

  float flop = 0.;

  {
    [[maybe_unused]] Profiler prf("Frequency FT: HOST", "tp-accumulation", __LINE__, thread_id_);
    for (int s = 0; s < 2; ++s)
      flop += ndft_objs_[stream_id(s)].execute(configs[s], M_pair[s], G_[s]);
  }
  {
    [[maybe_unused]] Profiler prf("Space FT: HOST", "tp-accumulation", __LINE__, thread_id_);
    for (int s = 0; s < 2; ++s)
      flop += space_trsf_objs_[stream_id(s)].execute(G_[s]);
  }

  return flop;
}

template <class Parameters, DistType DT>
void TpAccumulatorGpuBase<Parameters, DT>::sumTo_(TpAccumulatorGpuBase<Parameters, DT>& /*other_one*/) {
  // Nothing to do: G4 on the device is shared.
  synchronizeStreams();
  return;
}

}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif
