// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This class implements the 1D delayed-NFFT (d-NFFT) algorithm on the GPU.

#ifndef DCA_MATH_NFFT_DNFFT_1D_GPU_HPP
#define DCA_MATH_NFFT_DNFFT_1D_GPU_HPP
#ifdef DCA_HAVE_CUDA

#include "dca/math/nfft/dnfft_1d.hpp"

#include <cuda_runtime.h>

#include "dca/linalg/matrix.hpp"
#include "dca/linalg/util/allocators.hpp"
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

namespace details {
// dca::math::nfft::details::
using BDmn = func::dmn_0<phys::domains::electron_band_domain>;
template <class RDmn>
using PDmn = func::dmn_variadic<BDmn, BDmn, RDmn>;
}

// Only the CUBIC mode is implemented.
template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
class Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>
    : public Dnfft1D<ScalarType, WDmn, details::PDmn<RDmn>, oversampling, CUBIC> {
public:
  using BDmn = details::BDmn;
  using PDmn = details::PDmn<RDmn>;

  using ThisType = Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling>;
  using ElementType = ScalarType;
  using BaseClass = Dnfft1D<ScalarType, WDmn, PDmn, oversampling, CUBIC>;

public:
  Dnfft1DGpu(double beta, cudaStream_t stream, bool accumulate_m_sqr = false);

  void initialize();

  template <typename InpScalar, class Configuration>
  void accumulate(const linalg::Matrix<InpScalar, linalg::CPU>& M, const Configuration& config,
                  const int sign);

  template <typename OtherScalarType>
  void finalize(func::function<std::complex<OtherScalarType>, func::dmn_variadic<WDmn, PDmn>>& f_w,
                bool get_square = false);

  // Sums the accumulated data in the time domain.
  ThisType& operator+=(ThisType& other);

private:
  void initializeDeviceCoefficients();

  void uploadMatrix(const linalg::Matrix<ScalarType, linalg::CPU>& M);
  template <typename InpScalar>
  void uploadMatrix(const linalg::Matrix<InpScalar, linalg::CPU>& M);

private:
  using BaseClass::f_tau_;
  static linalg::Vector<ScalarType, linalg::GPU> cubic_coeff_dev_;

  const double beta_;
  cudaStream_t stream_;
  const bool accumulate_m_sqr_;
  linalg::Matrix<ScalarType, linalg::GPU> accumulation_matrix_;
  linalg::Matrix<ScalarType, linalg::GPU> accumulation_matrix_sqr_;

  linalg::Matrix<ScalarType, linalg::CPU> M_host_;
  linalg::Matrix<ScalarType, linalg::GPU> M_;
  linalg::Matrix<ScalarType, linalg::GPU> M_sqr_;
  linalg::util::HostVector<details::ConfigElem> config_left_;
  linalg::util::HostVector<details::ConfigElem> config_right_;
  linalg::util::HostVector<ScalarType> times_;
  linalg::Vector<details::ConfigElem, linalg::GPU> config_left_dev_;
  linalg::Vector<details::ConfigElem, linalg::GPU> config_right_dev_;
  linalg::Vector<ScalarType, linalg::GPU> times_dev_;

  linalg::util::CudaEvent config_copied_event_;
  linalg::util::CudaEvent m_copied_event_;
};
template <typename ScalarType, typename WDmn, typename PDmn, int oversampling>
linalg::Vector<ScalarType, linalg::GPU>
    Dnfft1DGpu<ScalarType, WDmn, PDmn, oversampling, CUBIC>::cubic_coeff_dev_;

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::Dnfft1DGpu(const double beta,
                                                                    cudaStream_t stream,
                                                                    const bool accumulate_m_sqr)
    : BaseClass(), beta_(beta), stream_(stream), accumulate_m_sqr_(accumulate_m_sqr) {
  initializeDeviceCoefficients();
  assert(cudaPeekAtLastError() == cudaSuccess);
  func::dmn_variadic<BDmn, BDmn, RDmn> bbr_dmn;
  const int n_times = BaseClass::PaddedTimeDmn::dmn_size();
  accumulation_matrix_.resizeNoCopy(std::make_pair(n_times, bbr_dmn.get_size()));
  if (accumulate_m_sqr)
    accumulation_matrix_sqr_.resizeNoCopy(accumulation_matrix_.size());
}

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
void Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::initialize() {
  accumulation_matrix_.setToZero(stream_);
  if (accumulate_m_sqr_)
    accumulation_matrix_sqr_.setToZero(stream_);
  assert(cudaPeekAtLastError() == cudaSuccess);
}

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
void Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::initializeDeviceCoefficients() {
  static bool initialized = false;
  static std::mutex mutex;
  if (initialized)
    return;

  std::unique_lock<std::mutex> lock(mutex);
  if (initialized)
    return;
  const auto& host_coeff = BaseClass::get_cubic_convolution_matrices();
  cubic_coeff_dev_.resizeNoCopy(host_coeff.size());
  cudaMemcpy(cubic_coeff_dev_.ptr(), host_coeff.values(), host_coeff.size() * sizeof(ScalarType),
             cudaMemcpyHostToDevice);

  const auto& sub_matrix = RDmn::parameter_type::get_subtract_matrix();
  using PaddedTimeDmn = typename BaseClass::PaddedTimeDmn::parameter_type;
  using WindowTimeDmn = typename BaseClass::WindowFunctionTimeDmn::parameter_type;
  details::initializeNfftHelper<ScalarType>(
      BDmn::dmn_size(), RDmn::dmn_size(), sub_matrix.ptr(), sub_matrix.leadingDimension(),
      oversampling, BaseClass::get_window_sampling(), PaddedTimeDmn::first_element(),
      PaddedTimeDmn::get_Delta(), WindowTimeDmn::first_element(), WindowTimeDmn::get_delta(), beta_);

  assert(cudaPeekAtLastError() == cudaSuccess);
  initialized = true;
}

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
template <typename InpScalar, class Configuration>
void Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::accumulate(
    const linalg::Matrix<InpScalar, linalg::CPU>& M, const Configuration& config, const int sign) {
  assert(M.is_square());
  const int n = M.nrCols();

  uploadMatrix(M);

  config_copied_event_.block();
  config_right_.resize(n);
  config_left_.resize(n);
  times_.resize(n);

  for (int i = 0; i < n; ++i) {
    // invert left and right as M is an inverse matrix.
    config_right_[i].b = config[i].get_left_band();
    config_right_[i].r = config[i].get_left_site();
    config_left_[i].b = config[i].get_right_band();
    config_left_[i].r = config[i].get_right_site();
    times_[i] = config[i].get_tau();
  }

  config_right_dev_.setAsync(config_right_, stream_);
  config_left_dev_.setAsync(config_left_, stream_);
  times_dev_.setAsync(times_, stream_);
  config_copied_event_.record(stream_);

  details::accumulateOnDevice(M_.ptr(), M_.leadingDimension(), sign, accumulation_matrix_.ptr(),
                              accumulation_matrix_sqr_.ptr(), accumulation_matrix_.leadingDimension(),
                              config_left_dev_.ptr(), config_right_dev_.ptr(), times_dev_.ptr(),
                              cubic_coeff_dev_.ptr(), n, stream_);
}

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
template <typename InpScalar>
void Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::uploadMatrix(
    const linalg::Matrix<InpScalar, linalg::CPU>& M) {
  m_copied_event_.block();
  M_host_.resizeNoCopy(M.size());
  for (int j = 0; j < M.nrCols(); ++j)
    for (int i = 0; i < M.nrRows(); ++i)
      M_host_(i, j) = static_cast<ScalarType>(M(i, j));
  M_.setAsync(M_host_, stream_);
  m_copied_event_.record(stream_);
}

template <typename ScalarType, typename WDmn, typename RDmn, int oversampling>
void Dnfft1DGpu<ScalarType, WDmn, RDmn, oversampling, CUBIC>::uploadMatrix(
    const linalg::Matrix<ScalarType, linalg::CPU>& M) {
  M_.setAsync(M, stream_);
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

}  // nfft
}  // math
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_MATH_NFFT_DNFFT_1D_GPU_HPP
