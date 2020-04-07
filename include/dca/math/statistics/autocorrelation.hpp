// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a class to compute the autocorrelation of a time series.

#ifndef DCA_MATH_STATISTICS_AUTOCORRELATION_HPP
#define DCA_MATH_STATISTICS_AUTOCORRELATION_HPP

#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dca/phys/dca_step/cluster_solver/shared_tools/util/accumulator.hpp"

namespace dca {
namespace math {
namespace statistics {

template <class T>
T conjugate(T x) {
  return x;
}
template <class T>
std::complex<T> conjugate(std::complex<T> x) {
  return std::conj(x);
}

template <class T>
class Autocorrelation {
public:
  using MeanType = typename phys::solver::util::Accumulator<T>::MeanType;

  // Constructs an object able to compute correlations up to a time difference of 't_max'.
  Autocorrelation(unsigned t_max = 0);

  void resize(std::size_t n);

  // Add a sample to the time-series, accumulating the autocorrelation on the fly.
  // Precondition: get() or computeAutocorrelationTime() have not been called from this object.
  void addSample(const T& sample);

  // Returns the autocorrelations at different times, with the first entry storing the variance.
  const auto& get();

  // Returns 1 + 2 \sum_{t >= 1} corr(t)
  auto computeAutocorrelationTime() -> MeanType;

  Autocorrelation& operator+=(const Autocorrelation& rhs);

  template <class Concurrency>
  void sumConcurrency(const Concurrency& concurrency);

  void reset();

  std::size_t samples() const {
    return mean_accum_.count();
  }

  auto getStdev() {
    return std::sqrt(get()[0]);
  }

  auto getMean() {
    return mean_accum_.mean();
  }

private:
  void finalize();

  bool finalized_ = false;
  std::vector<T> samples_;
  std::vector<MeanType> correlations_;
  std::vector<unsigned> counts_;
  unsigned samples_head_ = 0;
  phys::solver::util::Accumulator<T> mean_accum_;
};

template <class T>
Autocorrelation<T>::Autocorrelation(unsigned t_max) {
  resize(t_max);
}

template <class T>
void Autocorrelation<T>::resize(std::size_t n) {
  if (finalized_) {
    throw(std::logic_error("Object already finalized."));
  }
  samples_.resize(n);
  correlations_.resize(n, 0);
  counts_.resize(n, 0);
}

template <class T>
void Autocorrelation<T>::addSample(const T& sample) {
  if (finalized_) {
    throw(std::logic_error("Object already finalized."));
  }
  // Store sample
  samples_[samples_head_] = sample;
  mean_accum_.addSample(sample);

  ++samples_head_;
  if (samples_head_ == samples_.size()) {  // wrap around.
    samples_head_ = 0;
  }

  // accumulate correlations.
  unsigned n_accumulations = std::min(samples_.size(), mean_accum_.count());
  for (unsigned distance = 0; distance < n_accumulations; ++distance) {
    int idx = samples_head_ - 1 - distance;
    if (idx < 0)  // wrap-around
      idx += samples_.size();

    correlations_[distance] += samples_[idx] * conjugate(sample);
    counts_[distance]++;
  }
}

template <class T>
void Autocorrelation<T>::finalize() {
  if (finalized_)
    return;
  finalized_ = true;

  const auto mean = mean_accum_.mean();
  const auto mean2 = mean * conjugate(mean);

  correlations_[0] = correlations_[0] / static_cast<T>(mean_accum_.count()) - mean2;
  const auto c0 = correlations_[0];

  for (unsigned i = 1; i < correlations_.size(); ++i) {
    if (counts_[i] <= 0)
      break;
    correlations_[i] = (correlations_[i] / static_cast<T>(counts_[i]) - mean2) / c0;
  }
}

template <class T>
auto Autocorrelation<T>::computeAutocorrelationTime() -> MeanType {
  finalize();

  constexpr double factor = 6;

  bool converged = false;
  MeanType tau = 1.;

  for (unsigned t = 1; t < correlations_.size(); ++t) {
    if (counts_[t] <= 0)
      break;

    tau += static_cast<MeanType>(2. * (1. - double(t) / counts_[t])) * correlations_[t];

    if (t > factor * std::abs(tau)) {
      converged = true;
      break;
    }
  }

  if (converged == false) {
    std::cerr << "Autocorrelation time failed to converge." << std::endl;
  }

  return tau;
}

template <class T>
const auto& Autocorrelation<T>::get() {
  finalize();
  return correlations_;
}

template <class T>
Autocorrelation<T>& Autocorrelation<T>::operator+=(const Autocorrelation<T>& rhs) {
  if (correlations_.size() != rhs.correlations_.size()) {
    throw(std::logic_error("Autocorrelation size mismatch."));
  }

  mean_accum_ += rhs.mean_accum_;
  for (int i = 0; i < correlations_.size(); ++i) {
    correlations_[i] += rhs.correlations_[i];
    counts_[i] += rhs.counts_[i];
  }

  return *this;
}

template <class T>
template <class Concurrency>
void Autocorrelation<T>::sumConcurrency(const Concurrency& concurrency) {
  mean_accum_.sumConcurrency(concurrency);
  concurrency.sum(correlations_);
  concurrency.sum(counts_);
}

template <class T>
void Autocorrelation<T>::reset() {
  mean_accum_.reset();
  finalized_ = 0;
  samples_head_ = 0;
  std::fill(correlations_.begin(), correlations_.end(), 0);
  std::fill(samples_.begin(), samples_.end(), 0);
  std::fill(counts_.begin(), counts_.end(), 0);
}

}  // namespace statistics
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_STATISTICS_AUTOCORRELATION_HPP
