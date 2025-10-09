// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This class measures the single-particle functions with delayed NFFT scheme. The convolution is
// performed on the GPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_SP_SP_ACCUMULATOR_GPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_SP_SP_ACCUMULATOR_GPU_HPP

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/sp/sp_accumulator.hpp"

#include <cassert>
#include <complex>
#include <memory>
#include <stdexcept>
#include <vector>

#ifdef DCA_HAVE_GPU
#include "dca/platform/dca_gpu.h"
#include "dca/linalg/util/gpu_event.hpp"
#include "dca/linalg/util/gpu_stream.hpp"
#else
#error "This file requires GPU."
#endif

#include "dca/math/nfft/dnfft_1d_gpu.hpp"
#include "dca/phys/error_computation_type.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
// dca::phys::solver::accumulator::

template <class Parameters>
class SpAccumulator<Parameters, linalg::GPU>
    : public SpAccumulator<Parameters, linalg::CPU> {
public:
  using BaseClass = SpAccumulator<Parameters, linalg::CPU>;

  using typename BaseClass::BDmn;
  using typename BaseClass::PDmn;
  using typename BaseClass::RDmn;
  using typename BaseClass::WDmn;
  using typename BaseClass::Real;
  using typename BaseClass::Scalar;
  using typename BaseClass::Profiler;

  using BaseClass::accumulate_m_sqr_;
  using BaseClass::finalized_;
  using BaseClass::M_r_w_;
  using BaseClass::M_r_w_sqr_;
  using BaseClass::single_measurement_M_r_w_;
  using BaseClass::oversampling;
  using BaseClass::parameters_;
  using NfftType = math::nfft::Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling, math::nfft::CUBIC>;
  using MFunction = typename BaseClass::MFunction;
  using MFunctionTime = typename BaseClass::MFunctionTime;
  using MFunctionTimePair = typename BaseClass::MFunctionTimePair;
  using FTauPair = typename BaseClass::FTauPair;
  using FTau = typename BaseClass::FTau;
  using PaddedTimeDmn = typename NfftType::PaddedTimeDmn;
public:
  SpAccumulator(const Parameters& parameters_ref, bool accumulate_m_squared = false);

  void resetAccumulation();

  template <class Configuration>
  void accumulate(const std::array<linalg::Matrix<Scalar, linalg::GPU>, 2>& Ms,
                  const std::array<Configuration, 2>& configs, const Scalar factor);

  // For testing purposes.
  template <class Configuration>
  void accumulate(const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& Ms,
                  const std::array<Configuration, 2>& configs, const Scalar factor);

  void finalize();

  void sumTo(SpAccumulator<Parameters, linalg::GPU>& other);

  void synchronizeCopy() {
    cached_nfft_obj_[0].synchronizeCopy();
    cached_nfft_obj_[1].synchronizeCopy();
  }

  void syncStreams(const linalg::util::GpuEvent& event) {
    for (const auto& stream : streams_)
      event.block(stream);
  }

  auto get_streams() {
    return std::array<linalg::util::GpuStream*, 2>{&streams_[0], &streams_[0]};
  }

  // Returns the allocated device memory in bytes.
  int deviceFingerprint() const {
    return cached_nfft_obj_[0].deviceFingerprint() + cached_nfft_obj_[1].deviceFingerprint();
  }

  const MFunction& get_single_measurement_sign_times_MFunction();
  const FTauPair& get_single_measurement_sign_times_MFunction_time();

  template <class Writer>
  void write(Writer& writer) {
    writer.open_group(NfftType::PaddedTimeDmn::get_name());
    writer.execute("elements", NfftType::PaddedTimeDmn::get_elements());
    writer.close_group();
  }

private:
  /** finalize a std::array<NfftType, 2>
   *  \param[inout]      ft_objs
   *  \param[out]        function
   *  \param[in]         m_sqr      finalize the sqr elements of the NfftType(GPU)
   */
  void finalizeFunction(std::array<NfftType, 2>& ft_objs, MFunction& function, bool m_sqr);

  std::array<linalg::util::GpuStream, 1> streams_;
  /** gpu M_r_t */
  std::array<NfftType, 2> cached_nfft_obj_;
  /** \todo Don't always pay the memory cost even when not collect single measurement G's */
  std::array<NfftType, 2> single_measurement_M_r_t_device_;
};

/// \todo examine memory impact when not in use of single_measurement_M_r_w etc.

template <class Parameters>
SpAccumulator<Parameters, linalg::GPU>::SpAccumulator(const Parameters& parameters_ref,
                                                            const bool accumulate_m_sqr)
    : BaseClass(parameters_ref, accumulate_m_sqr),
      streams_(),
      cached_nfft_obj_{NfftType(parameters_.get_beta(), streams_[0], accumulate_m_sqr),
                       NfftType(parameters_.get_beta(), streams_[0], accumulate_m_sqr)},
      single_measurement_M_r_t_device_{NfftType(parameters_.get_beta(), streams_[0], false),
                                       NfftType(parameters_.get_beta(), streams_[0], false)} {
  single_measurement_M_r_w_.reset(new MFunction("M_r_w"));
}

template <class Parameters>
void SpAccumulator<Parameters, linalg::GPU>::resetAccumulation() {
  for (int s = 0; s < 2; ++s) {
    cached_nfft_obj_[s].resetAccumulation();
    if (parameters_.stamping_period() > 0) {
      single_measurement_M_r_t_device_[s].resetAccumulation();
      BaseClass::resetAccumulation();
      //single_measurement_M_r_t_->operator[](s).resetAccumulation();
    }
  }
  single_measurement_M_r_w_.reset(new MFunction("M_r_w"));
  finalized_ = false;
}

template <class Parameters>
template <class Configuration>
void SpAccumulator<Parameters, linalg::GPU>::accumulate(
    const std::array<linalg::Matrix<Scalar, linalg::GPU>, 2>& Ms,
    const std::array<Configuration, 2>& configs, const Scalar factor) {
  if (finalized_)
    throw(std::logic_error("The accumulator is already finalized."));

  for (int s = 0; s < 2; ++s) {
    cached_nfft_obj_[s].reserve(configs[s].size());
    if (parameters_.stamping_period() > 0) {
      // This is where I think I'm zeroing this out each time.
      single_measurement_M_r_t_device_[s].resetAccumulation();
      single_measurement_M_r_t_device_[s].reserve(configs[s].size());
    }
  }

  for (int s = 0; s < 2; ++s) {
    cached_nfft_obj_[s].accumulate(Ms[s], configs[s], factor);
    if (parameters_.stamping_period() > 0) {
      single_measurement_M_r_t_device_[s].accumulate(Ms[s], configs[s], factor);
    }
  }
}
  
template <class Parameters>
template <class Configuration>
void SpAccumulator<Parameters, linalg::GPU>::accumulate(
    const std::array<linalg::Matrix<Scalar, linalg::CPU>, 2>& Ms,
    const std::array<Configuration, 2>& configs, const Scalar factor) {
  std::array<linalg::Matrix<Scalar, linalg::GPU>, 2> M_dev;
  for (int s = 0; s < 2; ++s)
    M_dev[s].setAsync(Ms[s], streams_[0]);

  accumulate(M_dev, configs, factor);
}

template <class Parameters>
void SpAccumulator<Parameters, linalg::GPU>::finalizeFunction(std::array<NfftType, 2>& ft_objs,
                                                                    MFunction& function, bool m_sqr) {
  func::function<std::complex<Real>, func::dmn_variadic<WDmn, PDmn>> tmp("tmp");
  const Real normalization = 1. / RDmn::dmn_size();

  for (int s = 0; s < 2; ++s) {
    // This is the difference from CPU m_sqr falg because the sqr is bundled in nfft GPU object
    ft_objs[s].finalize(tmp, m_sqr);
    for (int w_ind = 0; w_ind < WDmn::dmn_size(); w_ind++)
      for (int r_ind = 0; r_ind < RDmn::dmn_size(); r_ind++)
        for (int b2_ind = 0; b2_ind < BDmn::dmn_size(); b2_ind++)
          for (int b1_ind = 0; b1_ind < BDmn::dmn_size(); b1_ind++)
            function(b1_ind, s, b2_ind, s, r_ind, w_ind) +=
                tmp(w_ind, b1_ind, b2_ind, r_ind) * normalization;
  }
}

template <class Parameters>
void SpAccumulator<Parameters, linalg::GPU>::finalize() {
  if (finalized_)
    return;

  M_r_w_.reset(new MFunction("M_r_w"));
  finalizeFunction(cached_nfft_obj_, *M_r_w_, false);

  if (accumulate_m_sqr_) {
    M_r_w_sqr_.reset(new MFunction("M_r_w_sqr"));
    finalizeFunction(cached_nfft_obj_, *M_r_w_sqr_, true);
  }

  finalized_ = true;
}

template <class Parameters>
void SpAccumulator<Parameters, linalg::GPU>::sumTo(
    SpAccumulator<Parameters, linalg::GPU>& other) {
  for (int s = 0; s < 2; ++s)
    other.cached_nfft_obj_[s] += cached_nfft_obj_[s];
}

// template <class Parameters>
// const typename SpAccumulator<Parameters, linalg::CPU>::MFunction& SpAccumulator<
//     Parameters, linalg::CPU>::get_single_measurement_sign_times_MFunction() {
//   single_measurement_M_r_w_.reset(new MFunction("single_function_M_r_w"));
//   finalizeFunction(*single_measurement_M_r_t_, *single_measurement_M_r_w_);
//   return *single_measurement_M_r_w_;
// }

/** get M_r_w_ for a single configuration.
 *  This is quite unoptimized, a heap allocation in incurred every time.
 */
template <class Parameters>
const typename SpAccumulator<Parameters, linalg::CPU>::MFunction& SpAccumulator<
    Parameters, linalg::GPU>::get_single_measurement_sign_times_MFunction() {
  // single_measurement_M_r_w_.reset(new MFunction("single_function_M_r_w"));
  //  assuming this is faster than the allocation.
  std::fill(single_measurement_M_r_w_->begin(), single_measurement_M_r_w_->end(),
            std::complex<double>{0.0, 0.0});
  finalizeFunction(single_measurement_M_r_t_device_, *single_measurement_M_r_w_, false);
  return *single_measurement_M_r_w_;
}

/** get M_r_t_ for a single configuration.
 *  This is quite unoptimized, a heap allocation in incurred every time.
 */
template <class Parameters>
const typename SpAccumulator<Parameters, linalg::CPU>::FTauPair& SpAccumulator<
    Parameters, linalg::GPU>::get_single_measurement_sign_times_MFunction_time() {
  // single_measurement_M_r_w_.reset(new MFunction("single_function_M_r_w"));
  //  assuming this is faster than the allocation.
  // std::fill(single_measurement_M_r_t_->begin(), single_measurement_M_r_t_->end(), 0.0);
  BaseClass::single_meas_ftau_pair_[0] = single_measurement_M_r_t_device_[0].get_f_tau();
  BaseClass::single_meas_ftau_pair_[1] = single_measurement_M_r_t_device_[1].get_f_tau();
  return BaseClass::single_meas_ftau_pair_;
}

}  // namespace accumulator
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_SP_SP_ACCUMULATOR_GPU_HPP
