// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This class implements the 1D delayed-NFFT (d-NFFT) algorithm on the GPU.
// See dnfft_1d.hpp for references and CPU implementation.

#ifndef DCA_MATH_NFFT_DNFFT_1D_GPU_HPP
#define DCA_MATH_NFFT_DNFFT_1D_GPU_HPP

#include "dca/config/haves_defines.hpp"
#ifdef DCA_HAVE_GPU
#include "dca/platform/dca_gpu.h"
#else
#error "This file requires GPU."
#endif

#include "dca/math/nfft/dnfft_1d.hpp"
#include "dca/linalg/util/copy.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/util/allocators/vectors_typedefs.hpp"
#include "dca/linalg/util/gpu_event.hpp"
#include "dca/linalg/vector.hpp"
#include "dca/math/nfft/kernels_interface.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

template <typename Scalar, typename WDmn, typename RDmn, int oversampling = 8, NfftModeNames mode = CUBIC>
class Dnfft1DGpu;

// Only the CUBIC mode is implemented.
template <typename Scalar, typename WDmn, typename RDmn, int oversampling>
class Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling, CUBIC>
    : public Dnfft1D<Scalar, WDmn,
                     func::dmn_variadic<func::dmn_0<phys::domains::electron_band_domain>,
                                        func::dmn_0<phys::domains::electron_band_domain>, RDmn>,
                     oversampling, CUBIC> {
public:
  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using PDmn = func::dmn_variadic<BDmn, BDmn, RDmn>;

  using ThisType = Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling>;
  using ElementType = Scalar;
  using BaseClass = Dnfft1D<Scalar, WDmn, PDmn, oversampling, CUBIC>;
  using FTau = typename BaseClass::FTau;
  using PaddedTimeDmn = typename BaseClass::PaddedTimeDmn;
  using Real = dca::util::RealAlias<Scalar>;
  Dnfft1DGpu(double beta, const linalg::util::GpuStream& stream, bool accumulate_m_sqr = false);

  // Resets the accumulated quantities. To be called before each DCA iteration.
  void resetAccumulation();

  // Accumulates asynchronously on the device the entries of the M matrix at times and orbitals
  // described by the provided configuration.
  // Postcondition: M and config shall not be modified until 'synchronizeCopy' is called.
  template <class Configuration>
  void accumulate(const linalg::Matrix<Scalar, linalg::GPU>& M, const Configuration& config,
                  const Scalar factor);

  // Transforms the accumulated data to the frequency domain and stores it in f_w.
  // If get_square is true, transforms the accumulated squared data instead.
  template <typename OtherScalar>
  void finalize(func::function<std::complex<OtherScalar>, func::dmn_variadic<WDmn, PDmn>>& f_w,
                bool get_square = false);

  // Preallocates memory for the configuration to avoid a crash during the execution of accumulate.
  void reserve(std::size_t size);

  // Sums the accumulated data in the time domain.
  ThisType& operator+=(ThisType& other);

  // Ensures that the arguments of the method 'accumulate' have been copied to the GPU.
  void synchronizeCopy() {
    m_copied_event_.block();
  }

  // Returns the allocated device memory in bytes.
  std::size_t deviceFingerprint() const {
    return accumulation_matrix_.deviceFingerprint() + accumulation_matrix_sqr_.deviceFingerprint();
  }

  static std::size_t staticDeviceFingerprint() {
    return get_device_cubic_coeff().deviceFingerprint();
  }

  FTau& get_f_tau() {
    getDeviceData(accumulation_matrix_);
    return BaseClass::get_f_tau();
  }

  const FTau& get_f_tau() const {
    getDeviceData(accumulation_matrix_);
    return BaseClass::get_f_tau();
  }

private:
  void initializeDeviceCoefficients();

  void getDeviceData(const linalg::Matrix<Scalar, linalg::GPU>& data);

  using BaseClass::f_tau_;
  static inline auto& get_device_cubic_coeff();

  const Real beta_;
  const linalg::util::GpuStream& stream_;
  const bool accumulate_m_sqr_;
  linalg::Matrix<Scalar, linalg::GPU> accumulation_matrix_;
  linalg::Matrix<Scalar, linalg::GPU> accumulation_matrix_sqr_;
  // Just what is this doing?
  linalg::Matrix<Scalar, linalg::CPU> accumulation_matrix_host_;

  linalg::util::HostVector<details::ConfigElem> config_left_;
  linalg::util::HostVector<details::ConfigElem> config_right_;
  linalg::util::HostVector<Real> times_;
  linalg::Vector<details::ConfigElem, linalg::GPU> config_left_dev_;
  linalg::Vector<details::ConfigElem, linalg::GPU> config_right_dev_;
  linalg::Vector<Real, linalg::GPU> times_dev_;

  linalg::util::GpuEvent m_copied_event_;
};

template <typename Scalar, typename WDmn, typename RDmn, int oversampling>
Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling, CUBIC>::Dnfft1DGpu(const double beta,
                                                                const linalg::util::GpuStream& stream,
                                                                const bool accumulate_m_sqr)
    : BaseClass(), beta_(beta), stream_(stream), accumulate_m_sqr_(accumulate_m_sqr) {
  initializeDeviceCoefficients();
  assert(cudaPeekAtLastError() == cudaSuccess);
  func::dmn_variadic<BDmn, BDmn, RDmn> bbr_dmn;
  const int n_times = BaseClass::PaddedTimeDmn::dmn_size();
  accumulation_matrix_.resizeNoCopy(std::make_pair(n_times, bbr_dmn.get_size()));

  // \todo make this optional it doubles size of host memeory
  accumulation_matrix_host_.resizeNoCopy(std::make_pair(n_times, bbr_dmn.get_size()));

  if (accumulate_m_sqr)
    accumulation_matrix_sqr_.resizeNoCopy(accumulation_matrix_.size());

  resetAccumulation();
}

template <typename Scalar, typename WDmn, typename RDmn, int oversampling>
void Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling, CUBIC>::getDeviceData(
    const linalg::Matrix<Scalar, linalg::GPU>& data) {
  checkRC(cudaMemcpy2DAsync(f_tau_.values(), f_tau_[0] * sizeof(Scalar),
                    data.ptr(), data.leadingDimension() * sizeof(Scalar),
                    data.nrRows() * sizeof(Scalar), data.nrCols(), cudaMemcpyDeviceToHost,
			    stream_));
  checkRC(cudaStreamSynchronize(stream_.streamActually()));
}

template <typename Scalar, typename WDmn, typename RDmn, int oversampling>
void Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling, CUBIC>::resetAccumulation() {
  accumulation_matrix_.setToZero(stream_);
  if (accumulate_m_sqr_)
    accumulation_matrix_sqr_.setToZero(stream_);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <typename Scalar, typename WDmn, typename RDmn, int oversampling>
void Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling, CUBIC>::initializeDeviceCoefficients() {
  static std::once_flag flag;
  std::call_once(flag, [&]() {
    const auto& host_coeff = BaseClass::get_cubic_convolution_matrices();
    auto& dev_coeff = get_device_cubic_coeff();
    dev_coeff.resizeNoCopy(host_coeff.size());
    checkRC(cudaMemcpy(dev_coeff.ptr(), host_coeff.values(), host_coeff.size() * sizeof(Real),
		       cudaMemcpyHostToDevice));

    const auto& sub_matrix = RDmn::parameter_type::get_subtract_matrix();
    const auto& add_matrix = RDmn::parameter_type::get_add_matrix();
    using PaddedTimeDmn = typename BaseClass::PaddedTimeDmn::parameter_type;
    using WindowTimeDmn = typename BaseClass::WindowFunctionTimeDmn::parameter_type;
    details::initializeNfftHelper(BDmn::dmn_size(), RDmn::dmn_size(), add_matrix.ptr(),
                                  add_matrix.leadingDimension(), sub_matrix.ptr(),
                                  sub_matrix.leadingDimension(), PaddedTimeDmn::first_element(),
                                  PaddedTimeDmn::get_Delta(), WindowTimeDmn::first_element(),
                                  WindowTimeDmn::get_delta(), beta_);

    assert(cudaPeekAtLastError() == cudaSuccess);
  });
}

template <typename Scalar, typename WDmn, typename RDmn, int oversampling>
void Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling, CUBIC>::reserve(std::size_t size) {
  config_right_.resize(size);
  config_left_.resize(size);
  times_.resize(size);
}

template <typename Scalar, typename WDmn, typename RDmn, int oversampling>
template <class Configuration>
void Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling, CUBIC>::accumulate(
    const linalg::Matrix<Scalar, linalg::GPU>& M, const Configuration& config, const Scalar factor) {
  assert(M.is_square());
  if (config.size() == 0)  // Contribution is zero.
    return;

  const int n = M.nrCols();
  config_right_.resize(n);
  config_left_.resize(n);
  times_.resize(n);

  for (int i = 0; i < n; ++i) {
    // invert left and right as M is an inverse matrix.
    config_right_[i].band = config[i].get_left_band();
    config_right_[i].site = config[i].get_left_site();
    config_left_[i].band = config[i].get_right_band();
    config_left_[i].site = config[i].get_right_site();
    times_[i] = config[i].get_tau();
  }

  config_right_dev_.setAsync(config_right_, stream_);
  config_left_dev_.setAsync(config_left_, stream_);
  times_dev_.setAsync(times_, stream_);

  details::accumulateOnDevice<oversampling, BaseClass::window_sampling_, Scalar, Real>(
      M.ptr(), M.leadingDimension(), factor, accumulation_matrix_.ptr(),
      accumulation_matrix_sqr_.ptr(), accumulation_matrix_.leadingDimension(), config_left_dev_.ptr(),
      config_right_dev_.ptr(), times_dev_.ptr(), get_device_cubic_coeff().ptr(), n, stream_);

  m_copied_event_.record(stream_.streamActually());
}

template <typename Scalar, typename WDmn, typename RDmn, int oversampling>
template <typename OtherScalar>
void Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling, CUBIC>::finalize(
    func::function<std::complex<OtherScalar>, func::dmn_variadic<WDmn, PDmn>>& f_w, bool get_square) {
  auto get_device_data = [&](const linalg::Matrix<Scalar, linalg::GPU>& data) {
    checkRC(cudaMemcpy2DAsync(f_tau_.values(), f_tau_[0] * sizeof(Scalar), data.ptr(),
                      data.leadingDimension() * sizeof(Scalar), data.nrRows() * sizeof(Scalar),
			      data.nrCols(), cudaMemcpyDeviceToHost, stream_));
    checkRC(cudaStreamSynchronize(stream_));
  };

  if (!get_square)
    get_device_data(accumulation_matrix_);
  else
    get_device_data(accumulation_matrix_sqr_);

  BaseClass::finalize(f_w);
}

template <typename Scalar, typename WDmn, typename RDmn, int oversampling>
Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling, CUBIC>& Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling,
                                                                CUBIC>::operator+=(ThisType& other) {
  checkRC(cudaStreamSynchronize(other.stream_));

  details::sum(other.accumulation_matrix_.ptr(), other.accumulation_matrix_.leadingDimension(),
               accumulation_matrix_.ptr(), accumulation_matrix_.leadingDimension(),
               accumulation_matrix_.nrRows(), accumulation_matrix_.nrCols(), stream_);
  if (accumulate_m_sqr_)
    details::sum(other.accumulation_matrix_sqr_.ptr(),
                 other.accumulation_matrix_sqr_.leadingDimension(), accumulation_matrix_sqr_.ptr(),
                 accumulation_matrix_sqr_.leadingDimension(), accumulation_matrix_sqr_.nrRows(),
                 accumulation_matrix_sqr_.nrCols(), stream_);

  checkRC(cudaStreamSynchronize(stream_));
  return *this;
}

template <typename Scalar, typename WDmn, typename RDmn, int oversampling>
auto& Dnfft1DGpu<Scalar, WDmn, RDmn, oversampling, CUBIC>::get_device_cubic_coeff() {
  static linalg::Vector<Real, linalg::GPU> coefficients;
  return coefficients;
}

}  // namespace nfft
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_NFFT_DNFFT_1D_GPU_HPP
