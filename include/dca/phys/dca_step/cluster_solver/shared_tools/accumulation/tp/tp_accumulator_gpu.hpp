// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Implementation of the two particle Green's function computation on the GPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_GPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_GPU_HPP
#ifdef DCA_HAVE_CUDA

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator.hpp"

#include <cuda.h>
#include <mutex>

#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/util/magma_queue.hpp"
#include "dca/math/function_transform/special_transforms/space_transform_2D_gpu.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/kernels_interface.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft_gpu.hpp"

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
  // In: accumulate_on_device: if false perform the last multiplication and accumulation step on the
  //                           host.
  // In: thread_id: thread id, only used by the profiler.
  TpAccumulator(
      const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
      const Parameters& pars, bool accumulate_on_device = true, int thread_id = 0);

  // Resets the object between DCA iterations.
  void initialize();

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

  void initializeG4Helpers();

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
  using BaseClass::q_idx_;
  using BaseClass::w_idx_;
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

  std::array<MatrixDev, 2> G0_;
  std::array<MatrixHost, 2> G0_host_;
  std::array<MatrixDev, 2> G_;
  std::array<MatrixHost, 2> G_matrix_host_;
  linalg::Vector<Complex, linalg::GPU> G4_dev_;

  bool accumulate_on_device_;
  bool finalized_ = false;
};

template <class Parameters>
TpAccumulator<Parameters, linalg::GPU>::TpAccumulator(
    const func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KDmn, WDmn>>& G0,
    const Parameters& pars, const bool accumulate_on_device, const int thread_id)
    : BaseClass(G0, pars, thread_id),
      queues_(),
      streams_{queues_[0].getStream(), queues_[1].getStream()},
      ndft_objs_{NdftType(queues_[0]), NdftType(queues_[1])},
      space_trsf_objs_{DftType(n_pos_frqs_, queues_[0]), DftType(n_pos_frqs_, queues_[1])},
      accumulate_on_device_(accumulate_on_device) {
  initializeG4Helpers();
  initialize();
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::initialize() {
  const int sp_index_offset = (WDmn::dmn_size() - WTpExtDmn::dmn_size()) / 2;
  static func::dmn_variadic<BDmn, KDmn, WTpExtDmn> bkw_dmn;
  static typename BaseClass::TpDomain tp_dmn;

  for (int s = 0; s < 2; ++s) {
    G0_host_[s].resizeNoCopy(std::make_pair(bkw_dmn.get_size(), BDmn::dmn_size()));
    for (int w = 0; w < WTpExtDmn::dmn_size(); ++w)
      for (int k = 0; k < KDmn::dmn_size(); ++k)
        for (int b2 = 0; b2 < n_bands_; ++b2)
          for (int b1 = 0; b1 < n_bands_; ++b1)
            G0_host_[s](bkw_dmn(b1, k, w), b2) = (*G0_ptr_)(b1, s, b2, s, k, w + sp_index_offset);

    G0_[s].setAsync(G0_host_[s], streams_[s]);
  }
  if (accumulate_on_device_)  // Try to allocate the memory for G4 on the device.
    try {
      G4_dev_.resizeNoCopy(tp_dmn.get_size());
      cudaMemsetAsync(G4_dev_.ptr(), 0, G4_dev_.size() * sizeof(Complex), streams_[0]);
    }
    catch (...) {
      G4_dev_.clear();
      accumulate_on_device_ = false;
    }

  finalized_ = false;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::initializeG4Helpers() {
  static bool initialized = false;
  static std::mutex mutex;
  std::unique_lock<std::mutex> lock(mutex);

  if (initialized)
    return;
  const auto& add_mat = KDmn::parameter_type::get_add_matrix();
  const auto& sub_mat = KDmn::parameter_type::get_subtract_matrix();
  details::initializeG4Helpers(n_bands_, KDmn::dmn_size(), WTpPosDmn::dmn_size(), q_idx_, w_idx_,
                               add_mat.ptr(), add_mat.leadingDimension(), sub_mat.ptr(),
                               sub_mat.leadingDimension());
  assert(cudaPeekAtLastError() == cudaSuccess);
  initialized = true;
}

template <class Parameters>
template <class Configuration, class Scalar>
void TpAccumulator<Parameters, linalg::GPU>::accumulate(
    const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& M_pair,
    const std::array<Configuration, 2>& configs, const int sign) {
  Profiler profiler("accumulate", "tp-accumulation", __LINE__, thread_id_);

  if (finalized_)
    throw(std::logic_error("The accumulator was already finalized."));

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
      if (!accumulate_on_device_)
        G_matrix_host_[s].setAsync(G_[s], streams_[s]);
    }
  }

  event_.record(streams_[1]);
  event_.block(streams_[0]);

  if (!accumulate_on_device_) {
    {
      Profiler prf("Wait on device", "tp-accumulation", __LINE__, thread_id_);
      cudaStreamSynchronize(streams_[0]);
    }
    auto& G_func = BaseClass::G_;
    func::dmn_variadic<BDmn, KDmn, WTpExtDmn> bkw_dmn;
    for (int w2 = 0; w2 < WTpExtDmn::dmn_size(); ++w2)
      for (int w1 = 0; w1 < WTpExtPosDmn::dmn_size(); ++w1)
        for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
          for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1)
            for (int s = 0; s < 2; ++s)
              for (int b2 = 0; b2 < n_bands_; ++b2)
                for (int b1 = 0; b1 < n_bands_; ++b1)
                  G_func(b1, b2, s, k1, k2, w1, w2) =
                      G_matrix_host_[s](bkw_dmn(b1, k1, w1), bkw_dmn(b2, k2, w2));
  }
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::computeGSingleband(const int s) {
  details::computeGSingleband(G_[s].ptr(), G_[s].leadingDimension(), G0_[s].ptr(), KDmn::dmn_size(),
                              n_pos_frqs_, beta_, streams_[s]);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::computeGMultiband(const int s) {
  details::computeGMultiband(G_[s].ptr(), G_[s].leadingDimension(), G0_[s].ptr(),
                             G0_[s].leadingDimension(), n_bands_, KDmn::dmn_size(), n_pos_frqs_,
                             beta_, streams_[s]);
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

  if (!accumulate_on_device_) {  // Accumulate on Host.
    BaseClass::updateG4();
    return;
  }

  switch (mode_) {
    case PARTICLE_HOLE_MAGNETIC:
      details::updateG4PHMagnetic(G4_dev_.ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
                                  G_[1].leadingDimension(), n_bands_, KDmn::dmn_size(),
                                  WTpPosDmn::dmn_size(), sign_, streams_[0]);
      break;
    case PARTICLE_HOLE_CHARGE:
      details::updateG4PHCharge(G4_dev_.ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
                                G_[1].leadingDimension(), n_bands_, KDmn::dmn_size(),
                                WTpPosDmn::dmn_size(), sign_, streams_[0]);
      break;
    case PARTICLE_HOLE_TRANSVERSE:
      details::updateG4PHTransv(G4_dev_.ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
                                G_[1].leadingDimension(), n_bands_, KDmn::dmn_size(),
                                WTpPosDmn::dmn_size(), sign_, streams_[0]);
      break;
    case PARTICLE_PARTICLE_UP_DOWN:
      details::updateG4PPUpDown(G4_dev_.ptr(), G_[0].ptr(), G_[0].leadingDimension(), G_[1].ptr(),
                                G_[1].leadingDimension(), n_bands_, KDmn::dmn_size(),
                                WTpPosDmn::dmn_size(), sign_, streams_[0]);
      break;
    default:
      throw(std::logic_error("Mode non supported."));
  }
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::finalize() {
  if (finalized_)
    return;

  if (accumulate_on_device_) {
    G4_.reset(new TpGreenFunction("G4"));
    cudaMemcpyAsync(G4_->values(), G4_dev_.ptr(), G4_->size() * sizeof(Complex),
                    cudaMemcpyDeviceToHost, streams_[0]);
    synchronize();
    G4_dev_.clear();
  }
  finalized_ = true;
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::synchronize() {
  for (auto stream : streams_)
    cudaStreamSynchronize(stream);
}

template <class Parameters>
void TpAccumulator<Parameters, linalg::GPU>::sumTo(this_type& other_one) {
  finalize();
  other_one.finalize();

  BaseClass::sumTo(other_one);
}

}  // accumulator
}  // solver
}  // phys
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TP_ACCUMULATOR_GPU_HPP
