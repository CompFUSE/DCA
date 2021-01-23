// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This class implements the 1D delayed-NFFT (d-NFFT) algorithm on the GPU.
// See dnfft_1d.hpp for references and CPU implementation.

#ifndef DCA_HAVE_CUDA
#error "This file requires CUDA."
#endif

#ifndef DCA_MATH_NFFT_DNFFT_1D_GPU_HPP
#define DCA_MATH_NFFT_DNFFT_1D_GPU_HPP

#include "dca/math/nfft/dnfft_1d.hpp"

#include <cuda_runtime.h>

#include "dca/linalg/matrix.hpp"
#include "dca/linalg/util/allocators/vectors_typedefs.hpp"
#include "dca/linalg/util/cuda_event.hpp"
#include "dca/linalg/vector.hpp"
#include "dca/math/nfft/kernels_interface.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling = 8, NfftModeNames mode = CUBIC>
class Dnfft1DGpu;

// Only the CUBIC mode is implemented.
template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
class Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>
    : public Dnfft1D<ScalarType, WDmn,
                     func::dmn_variadic<func::dmn_0<phys::domains::electron_band_domain>,
                                        func::dmn_0<phys::domains::electron_band_domain>, RDmn>,
                     oversampling, CUBIC> {
public:
  using BDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using PDmn = func::dmn_variadic<BDmn, BDmn, RDmn>;

  using ThisType = Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling>;
  using ElementType = ScalarType;
  using BaseClass = Dnfft1D<ScalarType, WDmn, PDmn, oversampling, CUBIC>;

  Dnfft1DGpu(double beta, const linalg::util::CudaStream& stream, bool accumulate_m_sqr = false);

  // Resets the accumulated quantities. To be called before each DCA iteration.
  void resetAccumulation();

  // Accumulates asynchronously on the device the entries of the M matrix at times and orbitals
  // described by the provided configuration.
  // Postcondition: M and config shall not be modified until 'synchronizeCopy' is called.
  template <class Configuration, typename RealInp>
  void accumulate(const linalg::Matrix<RealInp, linalg::GPU>& M, const Configuration& config,
                  const int sign);

  // Transforms the accumulated data to the frequency domain and stores it in f_w.
  // If get_square is true, transforms the accumulated squared data instead.
  template <typename OtherScalarType>
  void finalize(func::function<std::complex<OtherScalarType>, func::dmn_variadic<WDmn, PDmn>>& f_w,
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

private:
  void initializeDeviceCoefficients();

  using BaseClass::f_tau_;
  static inline linalg::Vector<ScalarType, linalg::GPU>& get_device_cubic_coeff();

  const double beta_;
  const linalg::util::CudaStream& stream_;
  const bool accumulate_m_sqr_;
  linalg::Matrix<ScalarType, linalg::GPU> accumulation_matrix_;
  linalg::Matrix<ScalarType, linalg::GPU> accumulation_matrix_sqr_;

  linalg::util::HostVector<details::ConfigElem> config_left_;
  linalg::util::HostVector<details::ConfigElem> config_right_;
  linalg::util::HostVector<ScalarType> times_;
  linalg::Vector<details::ConfigElem, linalg::GPU> config_left_dev_;
  linalg::Vector<details::ConfigElem, linalg::GPU> config_right_dev_;
  linalg::Vector<ScalarType, linalg::GPU> times_dev_;

  linalg::util::CudaEvent m_copied_event_;
};

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::Dnfft1DGpu(
    const double beta, const linalg::util::CudaStream& stream, const bool accumulate_m_sqr)
    : BaseClass(), beta_(beta), stream_(stream), accumulate_m_sqr_(accumulate_m_sqr) {
  initializeDeviceCoefficients();
  assert(cudaPeekAtLastError() == cudaSuccess);
  func::dmn_variadic<BDmn, BDmn, RDmn> bbr_dmn;
  const int n_times = BaseClass::PaddedTimeDmn::dmn_size();
  accumulation_matrix_.resizeNoCopy(std::make_pair(n_times, bbr_dmn.get_size()));
  if (accumulate_m_sqr)
    accumulation_matrix_sqr_.resizeNoCopy(accumulation_matrix_.size());

  resetAccumulation();
}

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
void Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::resetAccumulation() {
  accumulation_matrix_.setToZero(stream_);
  if (accumulate_m_sqr_)
    accumulation_matrix_sqr_.setToZero(stream_);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
void Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::initializeDeviceCoefficients() {
  static std::once_flag flag;
  std::call_once(flag, [&]() {
    const auto& host_coeff = BaseClass::get_cubic_convolution_matrices();
    auto& dev_coeff = get_device_cubic_coeff();
    dev_coeff.resizeNoCopy(host_coeff.size());
    cudaMemcpy(dev_coeff.ptr(), host_coeff.values(), host_coeff.size() * sizeof(ScalarType),
               cudaMemcpyHostToDevice);

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

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
void Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::reserve(std::size_t size) {
  config_right_.resize(size);
  config_left_.resize(size);
  times_.resize(size);
}

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
template <class Configuration, typename RealInp>
void Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::accumulate(
    const linalg::Matrix<RealInp, linalg::GPU>& M, const Configuration& config, const int sign) {
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

  details::accumulateOnDevice<oversampling, BaseClass::window_sampling_, RealInp, ScalarType>(
      M.ptr(), M.leadingDimension(), static_cast<ScalarType>(sign), accumulation_matrix_.ptr(),
      accumulation_matrix_sqr_.ptr(), accumulation_matrix_.leadingDimension(), config_left_dev_.ptr(),
      config_right_dev_.ptr(), times_dev_.ptr(), get_device_cubic_coeff().ptr(), n, stream_);

  m_copied_event_.record(stream_);
}

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
template <typename OtherScalarType>
void Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::finalize(
    func::function<std::complex<OtherScalarType>, func::dmn_variadic<WDmn, PDmn>>& f_w,
    bool get_square) {
  auto get_device_data = [&](const linalg::Matrix<ScalarType, linalg::GPU>& data) {
    cudaMemcpy2DAsync(f_tau_.values(), f_tau_[0] * sizeof(ScalarType), data.ptr(),
                      data.leadingDimension() * sizeof(ScalarType),
                      data.nrRows() * sizeof(ScalarType), data.nrCols(), cudaMemcpyDeviceToHost,
                      stream_);
    cudaStreamSynchronize(stream_);
  };

  if (!get_square)
    get_device_data(accumulation_matrix_);
  else
    get_device_data(accumulation_matrix_sqr_);

  BaseClass::finalize(f_w);
}

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>& Dnfft1DGpu<
    ScalarType, WDmn, RDmn, oversampling, CUBIC>::operator+=(ThisType& other) {
  cudaStreamSynchronize(other.stream_);

  details::sum(other.accumulation_matrix_.ptr(), other.accumulation_matrix_.leadingDimension(),
               accumulation_matrix_.ptr(), accumulation_matrix_.leadingDimension(),
               accumulation_matrix_.nrRows(), accumulation_matrix_.nrCols(), stream_);
  if (accumulate_m_sqr_)
    details::sum(other.accumulation_matrix_sqr_.ptr(),
                 other.accumulation_matrix_sqr_.leadingDimension(), accumulation_matrix_sqr_.ptr(),
                 accumulation_matrix_sqr_.leadingDimension(), accumulation_matrix_sqr_.nrRows(),
                 accumulation_matrix_sqr_.nrCols(), stream_);

  cudaStreamSynchronize(stream_);
  return *this;
}

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
linalg::Vector<ScalarType, linalg::GPU>& Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling,
                                                    CUBIC>::get_device_cubic_coeff() {
  static linalg::Vector<ScalarType, linalg::GPU> coefficients;
  return coefficients;
}

}  // namespace nfft
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_NFFT_DNFFT_1D_GPU_HPP
