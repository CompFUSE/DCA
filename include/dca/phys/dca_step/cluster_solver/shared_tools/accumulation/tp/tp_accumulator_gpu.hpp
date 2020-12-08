// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Implementation of the two particle Green's function computation on the GPU.

#ifndef DCA_HAVE_CUDA
#error "This file requires CUDA."
#endif

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_GPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_GPU_HPP

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator.hpp"

#include <cuda.h>
#include <mutex>
#include <vector>

#include "dca/config/mc_options.hpp"
#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/reshapable_matrix.hpp"
#include "dca/linalg/util/allocators/managed_allocator.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/util/magma_queue.hpp"
#include "dca/math/function_transform/special_transforms/space_transform_2D_gpu.hpp"
#include "dca/parallel/util/call_once_per_loop.hpp"
#include "dca/parallel/util/get_workload.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/g4_helper.cuh"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/kernels_interface.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft_gpu.hpp"
#include "dca/util/integer_division.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

template <class Parameters>
class TpAccumulator<Parameters, linalg::GPU> : public TpAccumulator<Parameters, linalg::CPU> {
protected:
  using this_type = TpAccumulator<Parameters, linalg::GPU>;
  using BaseClass = TpAccumulator<Parameters, linalg::CPU>;

  using RDmn = typename BaseClass::RDmn;
  using KDmn = typename BaseClass::KDmn;
  using NuDmn = typename BaseClass::NuDmn;
  using WDmn = typename BaseClass::WDmn;

  using typename BaseClass::Real;

public:
  // Constructor:
  // In: G0: non interacting greens function.
  // In: pars: parameters object.
  // In: thread_id: thread id, only used by the profiler.
  TpAccumulator(
      const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
      const Parameters& pars, int thread_id = 0);

  // Resets the object between DCA iterations.
  void resetAccumulation(unsigned int dca_loop);

  // Computes the two particles Greens function from the M matrix and accumulates it internally.
  // In: M_array: stores the M matrix for each spin sector.
  // In: configs: stores the walker's configuration for each spin sector.
  // In: sign: sign of the configuration.
  // Returns: number of flop.
  template <class Configuration, typename RealIn>
  float accumulate(const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M,
                   const std::array<Configuration, 2>& configs, int sign);

  // CPU input. For testing purposes.
  template <class Configuration>
  float accumulate(const std::array<linalg::Matrix<double, linalg::CPU>, 2>& M,
                   const std::array<Configuration, 2>& configs, int sign);

  // Downloads the accumulation result to the host.
  void finalize();

  // Sums on the host the accumulated Green's function to the accumulated Green's function of
  // other_acc.
  void sumTo(this_type& other_acc);

  const linalg::util::CudaStream* get_stream() {
    return &queues_[0].getStream();
  }

  void synchronizeCopy() {
    ndft_objs_[0].synchronizeCopy();
    ndft_objs_[1].synchronizeCopy();
  }

  void syncStreams(const linalg::util::CudaEvent& event) {
    for (const auto& stream : queues_)
      event.block(stream);
  }

  std::size_t deviceFingerprint() const {
    std::size_t res(0);
    for (const auto& work : workspaces_)
      res += work->deviceFingerprint();
    for (int s = 0; s < 2; ++s)
      res += ndft_objs_[s].deviceFingerprint() + G_[s].deviceFingerprint() +
             space_trsf_objs_[s].deviceFingerprint();
    return res;
  }

  static std::size_t staticDeviceFingerprint() {
    std::size_t res = 0;

    for (const auto& G4_channel : get_G4())
      res += G4_channel.deviceFingerprint();

    res += get_G0()[0].deviceFingerprint() + get_G0()[1].deviceFingerprint();

    return res;
  }

protected:
  using Profiler = typename Parameters::profiler_type;

  using typename BaseClass::WTpDmn;
  using typename BaseClass::WTpExtDmn;
  using typename BaseClass::WTpExtPosDmn;
  using typename BaseClass::WTpPosDmn;

  using typename BaseClass::BDmn;
  using typename BaseClass::SDmn;

  using typename BaseClass::Complex;
  using typename BaseClass::TpGreensFunction;

  using Matrix = linalg::Matrix<Complex, linalg::GPU>;

  void initializeG0();
  void resetG4();
  void initializeG4Helpers() const;

  void computeG();

  void computeGMultiband(int s);

  void computeGSingleband(int s);

  template <class Configuration, typename RealIn>
  float computeM(const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M_pair,
                 const std::array<Configuration, 2>& configs);

  float updateG4(const std::size_t channel_index);

  void synchronizeStreams();

  constexpr static int n_ndft_queues_ = config::McOptions::memory_savings ? 1 : 2;

  using BaseClass::beta_;
  using BaseClass::channels_;
  using BaseClass::extension_index_offset_;
  using BaseClass::G0_ptr_;
  using BaseClass::G4_;
  using BaseClass::multiple_accumulators_;
  using BaseClass::n_bands_;
  using BaseClass::n_pos_frqs_;

  using BaseClass::non_density_density_;
  using BaseClass::sign_;
  using BaseClass::thread_id_;

  using MatrixDev = linalg::Matrix<Complex, linalg::GPU>;
  using RMatrix =
      linalg::ReshapableMatrix<Complex, linalg::GPU, config::McOptions::TpAllocator<Complex>>;
  using RMatrixValueType = typename RMatrix::ValueType;

  using MatrixHost = linalg::Matrix<Complex, linalg::CPU>;

  std::array<linalg::util::MagmaQueue, 2> queues_;
  linalg::util::CudaEvent event_;

  std::vector<std::shared_ptr<RMatrix>> workspaces_;

  using NdftType = CachedNdft<Real, RDmn, WTpExtDmn, WTpExtPosDmn, linalg::GPU, non_density_density_>;
  std::array<NdftType, 2> ndft_objs_;
  using DftType = math::transform::SpaceTransform2DGpu<RDmn, KDmn, Real>;
  std::array<DftType, 2> space_trsf_objs_;

  std::array<RMatrix, 2> G_;

  const int nr_accumulators_;

  bool finalized_ = false;
  bool initialized_ = false;

  using G0DevType = std::array<MatrixDev, 2>;
  static inline G0DevType& get_G0();
  using G4DevType = linalg::Vector<Complex, linalg::GPU, config::McOptions::TpAllocator<Complex>>;
  static inline std::vector<G4DevType>& get_G4();
};

template <class Parameters>
TpAccumulator<Parameters, linalg::GPU>::TpAccumulator(
    const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
    const Parameters& pars, const int thread_id)
    : BaseClass(G0, pars, thread_id),
      queues_(),
      ndft_objs_{NdftType(queues_[0]), NdftType(queues_[1])},
      space_trsf_objs_{DftType(n_pos_frqs_, queues_[0]), DftType(n_pos_frqs_, queues_[1])},
      nr_accumulators_(pars.get_accumulators()) {
  initializeG4Helpers();

  // Create shared workspaces.
  for (int i = 0; i < n_ndft_queues_; ++i) {
    workspaces_.emplace_back(std::make_shared<RMatrix>());
    workspaces_[i]->setStream(queues_[i]);
    ndft_objs_[i].setWorkspace(workspaces_[i]);
    space_trsf_objs_[i].setWorkspace(workspaces_[i]);
  }
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::resetAccumulation(const unsigned int dca_loop) {
  static dca::util::OncePerLoopFlag flag;

  dca::util::callOncePerLoop(flag, dca_loop, [&]() {
    resetG4();
    initializeG0();
    synchronizeStreams();
  });

  initialized_ = true;
  finalized_ = false;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::initializeG0() {
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

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::resetG4() {
  // Note: this method is not thread safe by itself.
  get_G4().resize(G4_.size());

  for (auto& G4_channel : get_G4()) {
    try {
      typename BaseClass::TpDomain tp_dmn;
      if (!multiple_accumulators_) {
        G4_channel.setStream(queues_[0].getStream());
      }
      G4_channel.resizeNoCopy(tp_dmn.get_size());
      G4_channel.setToZeroAsync(queues_[0].getStream());
    }
    catch (std::bad_alloc& err) {
      std::cerr << "Failed to allocate G4 on device.\n";
      throw(err);
    }
  }
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::initializeG4Helpers() const {
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

template <class Parameters>
template <class Configuration, typename RealIn>
float TpAccumulator<Parameters, linalg::GPU, DistType::NONE>::accumulate(
    const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M,
    const std::array<Configuration, 2>& configs, const int sign) {
  Profiler profiler("accumulate", "tp-accumulation", __LINE__, thread_id_);
  float flop = 0;

  if (!initialized_)
    throw(std::logic_error("The accumulator is not ready to measure."));

  if (!(configs[0].size() + configs[0].size()))  // empty config
    return flop;

  sign_ = sign;
  flop += computeM(M, configs);
  computeG();

  for (std::size_t channel = 0; channel < G4_.size(); ++channel) {
    flop += updateG4(channel);
  }

  return flop;
}

template <class Parameters>
template <class Configuration>
float TpAccumulator<Parameters, linalg::GPU>::accumulate(
    const std::array<linalg::Matrix<double, linalg::CPU>, 2>& M,
    const std::array<Configuration, 2>& configs, const int sign) {
  std::array<linalg::Matrix<double, linalg::GPU>, 2> M_dev;
  for (int s = 0; s < 2; ++s)
    M_dev[s].setAsync(M[s], queues_[0].getStream());

  return accumulate(M_dev, configs, sign);
}

template <class Parameters>
template <class Configuration, typename RealIn>
float TpAccumulator<Parameters, linalg::GPU>::computeM(
    const std::array<linalg::Matrix<RealIn, linalg::GPU>, 2>& M_pair,
    const std::array<Configuration, 2>& configs) {
  auto stream_id = [&](const int s) { return n_ndft_queues_ == 1 ? 0 : s; };

  float flop = 0.;

  {
    Profiler prf("Frequency FT: HOST", "tp-accumulation", __LINE__, thread_id_);
    for (int s = 0; s < 2; ++s)
      flop += ndft_objs_[stream_id(s)].execute(configs[s], M_pair[s], G_[s]);
  }
  {
    Profiler prf("Space FT: HOST", "tp-accumulation", __LINE__, thread_id_);
    for (int s = 0; s < 2; ++s)
      flop += space_trsf_objs_[stream_id(s)].execute(G_[s]);
  }

  return flop;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::computeG() {
  if (n_ndft_queues_ == 1) {
    event_.record(queues_[0]);
    event_.block(queues_[1]);
  }
  {
    Profiler prf("ComputeG: HOST", "tp-accumulation", __LINE__, thread_id_);
    for (int s = 0; s < 2; ++s) {
      if (n_bands_ == 1)
        computeGSingleband(s);
      else
        computeGMultiband(s);
    }
  }

  event_.record(queues_[1]);
  event_.block(queues_[0]);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::computeGSingleband(const int s) {
  details::computeGSingleband(G_[s].ptr(), G_[s].leadingDimension(), get_G0()[s].ptr(),
                              KDmn::dmn_size(), n_pos_frqs_, beta_, queues_[s]);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::computeGMultiband(const int s) {
  details::computeGMultiband(G_[s].ptr(), G_[s].leadingDimension(), get_G0()[s].ptr(),
                             get_G0()[s].leadingDimension(), n_bands_, KDmn::dmn_size(),
                             n_pos_frqs_, beta_, queues_[s]);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <class Parameters>
float TpAccumulator<Parameters, linalg::GPU>::updateG4(const std::size_t channel_index) {
  // G4 is stored with the following band convention:
  // b1 ------------------------ b3
  //        |           |
  //        |           |
  //        |           |
  // b2 ------------------------ b4

  //  TODO: set stream only if this thread gets exclusive access to G4.
  //  get_G4().setStream(queues_[0]);
  const FourPointType channel = channels_[channel_index];
  typename BaseClass::TpDomain tp_dmn;
  uint64_t start = 0;
  uint64_t end = tp_dmn.get_size();
  switch (channel) {
    case PARTICLE_HOLE_TRANSVERSE:
      return details::updateG4<Real, PARTICLE_HOLE_TRANSVERSE>(
          get_G4()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    case PARTICLE_HOLE_MAGNETIC:
      return details::updateG4<Real, PARTICLE_HOLE_MAGNETIC>(
          get_G4()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    case PARTICLE_HOLE_CHARGE:
      return details::updateG4<Real, PARTICLE_HOLE_CHARGE>(
          get_G4()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    case PARTICLE_HOLE_LONGITUDINAL_UP_UP:
      return details::updateG4<Real, PARTICLE_HOLE_LONGITUDINAL_UP_UP>(
          get_G4()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    case PARTICLE_HOLE_LONGITUDINAL_UP_DOWN:
      return details::updateG4<Real, PARTICLE_HOLE_LONGITUDINAL_UP_DOWN>(
          get_G4()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    case PARTICLE_PARTICLE_UP_DOWN:
      return details::updateG4<Real, PARTICLE_PARTICLE_UP_DOWN>(
          get_G4()[channel_index].ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), sign_, multiple_accumulators_, queues_[0], start, end);

    default:
      throw std::logic_error("Specified four point type not implemented.");
  }
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::finalize() {
  if (finalized_)
    return;

  for (std::size_t channel = 0; channel < G4_.size(); ++channel) {
    get_G4()[channel].copyTo(G4_[channel]);
  }
  // TODO: release memory if needed by the rest of the DCA loop.
  // get_G4().clear();

  finalized_ = true;
  initialized_ = false;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::synchronizeStreams() {
  for (auto& stream : queues_)
    cudaStreamSynchronize(stream);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::sumTo(this_type& /*other_one*/) {
  // Nothing to do: G4 on the device is shared.
  synchronizeStreams();
  return;
}

template <class Parameters>
auto TpAccumulator<Parameters, linalg::GPU>::get_G0() -> G0DevType& {
  static G0DevType G0;
  return G0;
}

template <class Parameters>
auto TpAccumulator<Parameters, linalg::GPU>::get_G4() -> std::vector<G4DevType>& {
  static std::vector<G4DevType> G4;
  return G4;
}

}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_GPU_HPP
