// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class implements a generic Monte Carlo accumulator.
//
// The template parameter T of Accumulator must either be
// - a floating point type (std::is_floating_point<T>::value == true),
// - an integer type (std::is_integral<T>::value == true),
// - or one of the floating point specializations of std::complex, i.e. std::complex<float>,
//   std::complex<double>, or std::complex<long double>.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_UTIL_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_UTIL_ACCUMULATOR_HPP

#include <cstdlib>  // std::size_t
#include <complex>
#include <stdexcept>
#include <type_traits>  // std::enable_if, std::is_floating_point, std::is_integral

namespace dca {
namespace phys {
namespace solver {
namespace util {
namespace details {
// dca::phys::solver::util::details::

// Primary template without member type.
template <typename T, typename Enable = void>
struct MeanType {};

// Specialization for floating point types.
template <typename T>
struct MeanType<T, typename std::enable_if_t<std::is_floating_point<T>::value>> {
  using type = T;
};

// Specialization for integer types.
template <typename T>
struct MeanType<T, typename std::enable_if_t<std::is_integral<T>::value>> {
  using type = double;
};

// Specialization for std::complex.
template <typename T>
struct MeanType<std::complex<T>, typename std::enable_if_t<std::is_floating_point<T>::value>> {
  using type = std::complex<T>;
};

}  // namespace details

template <typename T>
class Accumulator {
public:
  using SampleType = T;
  using MeanType = typename details::MeanType<T>::type;
  using CountType = std::size_t;

  Accumulator() : count_{}, sum_{} {}

  void addSample(const SampleType& sample) {
    sum_ += sample;
    ++count_;
  }

  Accumulator& operator+=(const Accumulator& other) {
    count_ += other.count_;
    sum_ += other.sum_;
    return *this;
  }

  CountType count() const {
    return count_;
  }

  SampleType sum() const {
    return sum_;
  }

  // Returns the current sample mean.
  // Preconditions: addSample has been called at least once after construction or last reset.
  MeanType mean() const {
    if (count_ == 0)
      throw std::logic_error("Empty sample has no mean.");

    return MeanType(sum_) / MeanType(count_);
  }

  void reset() {
    count_ = {};
    sum_ = {};
  }

  template <class Concurrency>
  void sumConcurrency(const Concurrency& concurrency) {
    // TODO: delay
    concurrency.sum(count_);
    concurrency.sum(sum_);
  }

private:
  CountType count_;
  SampleType sum_;
};

}  // namespace util
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_UTIL_ACCUMULATOR_HPP
