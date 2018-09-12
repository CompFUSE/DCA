// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
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

#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/util/magma_queue.hpp"
#include "dca/math/function_transform/special_transforms/space_transform_2D_gpu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/kernels_interface.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft_gpu.hpp"
#include "dca/util/call_once_per_loop.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

template <class Parameters>
class TpAccumulator<Parameters, linalg::GPU> : public TpAccumulator<Parameters, linalg::CPU> {
private:
  using this_type = TpAccumulator<Parameters, linalg::GPU>;
  using BaseClass = TpAccumulator<Parameters, linalg::CPU>;

  using RDmn = typename BaseClass::RDmn;
  using KDmn = typename BaseClass::KDmn;
  using NuDmn = typename BaseClass::NuDmn;
  using WDmn = typename BaseClass::WDmn;

public:
  // Constructor:
  // In: G0: non interacting greens function.
  // In: pars: parameters object.
  // In: thread_id: thread id, only used by the profiler.
  TpAccumulator(
      const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
      const Parameters& pars, int thread_id = 0);

  // Resets the object between DCA iterations.
  void resetAccumulation(uint dca_loop = 0);

  // Computes the two particles Greens function from the M matrix and accumulates it internally.
  // In: M_array: stores the M matrix for each spin sector.
  // In: configs: stores the walker's configuration for each spin sector.
  // In: sign: sign of the configuration.
  template <class Configuration, class Scalar>
  void accumulate(const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& M_pair,
                  const std::array<Configuration, 2>& configs, int sign);

  // Downloads the accumulation result to the host.
  void finalize();

  // Sums on the host the accumulated Green's function to the accumulated Green's function of
  // other_acc.
  void sumTo(this_type& other_acc);

  auto get_stream() const {
    return streams_[0];
  }

  void synchronizeCopy() {
    ndft_objs_[0].synchronizeCopy();
    ndft_objs_[1].synchronizeCopy();
  }

private:
  using Profiler = typename Parameters::profiler_type;

  using typename BaseClass::WTpDmn;
  using typename BaseClass::WTpPosDmn;
  using typename BaseClass::WTpExtDmn;
  using typename BaseClass::WTpExtPosDmn;

  using typename BaseClass::BDmn;
  using typename BaseClass::SDmn;

  using typename BaseClass::TpGreenFunction;
  using typename BaseClass::Real;
  using typename BaseClass::Complex;

  void initializeG0();
  void resetG4();
  void initializeG4Helpers() const;

  void computeG();

  void computeGMultiband(int s);

  void computeGSingleband(int s);

  template <class Configuration, typename Scalar>
  void computeM(const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& M_pair,
                const std::array<Configuration, 2>& configs);

  void updateG4();

  void synchronize();

private:
  using BaseClass::thread_id_;
  using BaseClass::n_bands_;
  using BaseClass::beta_;
  using BaseClass::G0_ptr_;
  using BaseClass::G4_;
  using BaseClass::mode_;
  using BaseClass::non_density_density_;
  using BaseClass::sign_;
  using BaseClass::extension_index_offset_;
  using BaseClass::n_pos_frqs_;

  using MatrixDev = linalg::Matrix<Complex, linalg::GPU>;
  using MatrixHost = linalg::Matrix<Complex, linalg::CPU>;

  std::array<linalg::util::MagmaQueue, 2> queues_;
  std::array<cudaStream_t, 2> streams_;
  linalg::util::CudaEvent event_;

  using NdftType = CachedNdft<Real, RDmn, WTpExtDmn, WTpExtPosDmn, linalg::GPU, non_density_density_>;
  std::array<NdftType, 2> ndft_objs_;
  using DftType = math::transform::SpaceTransform2DGpu<RDmn, KDmn, Real>;
  std::array<DftType, 2> space_trsf_objs_;

  std::array<MatrixDev, 2> G_;
  std::array<MatrixHost, 2> G_matrix_host_;

  bool finalized_ = false;
  bool initialized_ = false;

  using G0DevType = std::array<MatrixDev, 2>;
  static inline G0DevType& get_G0();
  using G4DevType = linalg::Vector<Complex, linalg::GPU>;
  static inline G4DevType& get_G4();
};

template <class Parameters>
TpAccumulator<Parameters, linalg::GPU>::TpAccumulator(
    const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
    const Parameters& pars, const int thread_id)
    : BaseClass(G0, pars, thread_id),
      queues_(),
      streams_{queues_[0].getStream(), queues_[1].getStream()},
      ndft_objs_{NdftType(queues_[0]), NdftType(queues_[1])},
      space_trsf_objs_{DftType(n_pos_frqs_, queues_[0]), DftType(n_pos_frqs_, queues_[1])} {
  initializeG4Helpers();
  resetAccumulation();
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::resetAccumulation(const uint dca_loop) {
  static util::OncePerLoopFlag flag;

  util::callOncePerLoop(flag, dca_loop, [&]() {
    resetG4();
    initializeG0();
    synchronize();
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

    G0[s].setAsync(G0_host[s], streams_[s]);
  }
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::resetG4() {
  // Note: this method is not thread safe by itself.
  auto& G4 = get_G4();
  try {
    typename BaseClass::TpDomain tp_dmn;
    G4.resizeNoCopy(tp_dmn.get_size());
    cudaMemsetAsync(G4.ptr(), 0, G4.size() * sizeof(Complex), streams_[0]);
  }
  catch (std::bad_alloc& err) {
    std::cerr << "Failed to allocate G4 on device.\n";
    throw(err);
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
    details::initializeG4Helpers(n_bands_, KDmn::dmn_size(), WTpPosDmn::dmn_size(), q_indices,
                                 w_indices, add_mat.ptr(), add_mat.leadingDimension(),
                                 sub_mat.ptr(), sub_mat.leadingDimension());
    assert(cudaPeekAtLastError() == cudaSuccess);
  });
}

template <class Parameters>
template <class Configuration, class Scalar>
void TpAccumulator<Parameters, linalg::GPU>::accumulate(
    const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& M_pair,
    const std::array<Configuration, 2>& configs, const int sign) {
  Profiler profiler("accumulate", "tp-accumulation", __LINE__, thread_id_);

  if (!initialized_)
    throw(std::logic_error("The accumulator is not ready to measure."));

  if (!(configs[0].size() + configs[0].size()))  // empty config
    return;

  sign_ = sign;
  computeM(M_pair, configs);
  computeG();
  updateG4();
}

template <class Parameters>
template <class Configuration, class Scalar>
void TpAccumulator<Parameters, linalg::GPU>::computeM(
    const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& M_pair,
    const std::array<Configuration, 2>& configs) {
  {
    Profiler prf("Frequency FT: HOST", "tp-accumulation", __LINE__, thread_id_);
    for (int s = 0; s < 2; ++s)
      ndft_objs_[s].execute(configs[s], M_pair[s], G_[s]);
  }
  {
    Profiler prf("Space FT: HOST", "tp-accumulation", __LINE__, thread_id_);
    for (int s = 0; s < 2; ++s)
      space_trsf_objs_[s].execute(G_[s]);
  }
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::computeG() {
  {
    Profiler prf("ComputeG: HOST", "tp-accumulation", __LINE__, thread_id_);
    for (int s = 0; s < 2; ++s) {
      if (n_bands_ == 1)
        computeGSingleband(s);
      else
        computeGMultiband(s);
    }
  }

  event_.record(streams_[1]);
  event_.block(streams_[0]);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::computeGSingleband(const int s) {
  details::computeGSingleband(G_[s].ptr(), G_[s].leadingDimension(), get_G0()[s].ptr(),
                              KDmn::dmn_size(), n_pos_frqs_, beta_, streams_[s]);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::computeGMultiband(const int s) {
  details::computeGMultiband(G_[s].ptr(), G_[s].leadingDimension(), get_G0()[s].ptr(),
                             get_G0()[s].leadingDimension(), n_bands_, KDmn::dmn_size(),
                             n_pos_frqs_, beta_, streams_[s]);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::updateG4() {
  // G4 is stored with the following band convention:
  // b1 ------------------------ b3
  //        |           |
  //        |           |
  //        |           |
  // b2 ------------------------ b4

  const int nw_exchange = domains::FrequencyExchangeDomain::get_size();
  const int nk_exchange = domains::MomentumExchangeDomain::get_size();

  switch (mode_) {
    case PARTICLE_HOLE_MAGNETIC:
      details::updateG4<Real, PARTICLE_HOLE_MAGNETIC>(
          get_G4().ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), n_bands_, KDmn::dmn_size(), WTpPosDmn::dmn_size(), nw_exchange,
          nk_exchange, sign_, streams_[0]);
      break;
    case PARTICLE_HOLE_CHARGE:
      details::updateG4<Real, PARTICLE_HOLE_CHARGE>(
          get_G4().ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), n_bands_, KDmn::dmn_size(), WTpPosDmn::dmn_size(), nw_exchange,
          nk_exchange, sign_, streams_[0]);
      break;
    case PARTICLE_HOLE_TRANSVERSE:
      details::updateG4<Real, PARTICLE_HOLE_TRANSVERSE>(
          get_G4().ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), n_bands_, KDmn::dmn_size(), WTpPosDmn::dmn_size(), nw_exchange,
          nk_exchange, sign_, streams_[0]);
      break;
    case PARTICLE_PARTICLE_UP_DOWN:
      details::updateG4<Real, PARTICLE_PARTICLE_UP_DOWN>(
          get_G4().ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
          G_[1].leadingDimension(), n_bands_, KDmn::dmn_size(), WTpPosDmn::dmn_size(), nw_exchange,
          nk_exchange, sign_, streams_[0]);
      break;
    default:
      throw(std::logic_error("Mode non supported."));
  }
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::finalize() {
  if (finalized_)
    return;

  G4_ = std::make_unique<TpGreenFunction>("G4");
  assert(G4_->size() == get_G4().size());

  cudaMemcpyAsync(G4_->values(), get_G4().ptr(), G4_->size() * sizeof(Complex),
                  cudaMemcpyDeviceToHost, streams_[0]);
  synchronize();
  // TODO: release memory if needed by the rest of the DCA loop.
  // get_G4().clear();

  finalized_ = true;
  initialized_ = false;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::synchronize() {
  for (auto stream : streams_)
    cudaStreamSynchronize(stream);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::sumTo(this_type& /*other_one*/) {
  // Nothing to do: G4 on the device is shared.
  synchronize();
  return;
}

template <class Parameters>
typename TpAccumulator<Parameters, linalg::GPU>::G0DevType& TpAccumulator<Parameters,
                                                                          linalg::GPU>::get_G0() {
  static G0DevType G0;
  return G0;
}

template <class Parameters>
typename TpAccumulator<Parameters, linalg::GPU>::G4DevType& TpAccumulator<Parameters,
                                                                          linalg::GPU>::get_G4() {
  static G4DevType G4;
  return G4;
}

}  // accumulator
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_GPU_HPP
