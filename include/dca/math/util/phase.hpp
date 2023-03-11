// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Class to represent and update either the sign of a real value, or the phase of a complex number.

#ifndef DCA_MATH_UTIL_PHASE_HPP
#define DCA_MATH_UTIL_PHASE_HPP

#include <cmath>
#include <complex>
#include <type_traits>

#include "dca/util/type_utils.hpp"

namespace dca::math {

template <bool is_complex>
class PhaseImpl;
  
template <>
class PhaseImpl<true> {
public:
  PhaseImpl() noexcept = default;

  template <class T>
  PhaseImpl(const std::complex<T>& val) noexcept {
    multiply(val);
  }

  template <class T>
  void multiply(const T& val) noexcept {
    if (val == T(0.))
      makeNull();
    phase_ += std::arg(val);
    adjustPhase();
  }

  void multiply(const PhaseImpl& other) noexcept {
    phase_ += other.phase_;
    adjustPhase();
  }

  template <class T>
  void divide(const T& val) noexcept {
    if (val == T(0.))
      makeNull();
    phase_ -= std::arg(val);
    adjustPhase();
  }

  void divide(const PhaseImpl& other) noexcept {
    phase_ -= other.phase_;
    adjustPhase();
  }

  void flip() noexcept {
    phase_ += M_PI;
    adjustPhase();
  }

  void reset() noexcept {
    phase_ = 0;
  }

  void makeNull() noexcept {
    phase_ = std::numeric_limits<double>::quiet_NaN();
  }

  bool isNull() const noexcept {
    return std::isnan(phase_);
  };

  std::complex<double> getSign() const noexcept {
    if (isNull())
      return {0, 0};
    else
      return std::polar(1., phase_);
  }

private:
  // map phase_ back into [-2*pi, 2*pi] (redundancy is intended).
  void adjustPhase() noexcept {
    phase_ = std::fmod(phase_, 2 * M_PI);
  }

  double phase_ = 0.;
};

template <>
class PhaseImpl<false> {
public:
  static constexpr bool IsPhase = true;
  PhaseImpl() noexcept = default;

  template <class T>
  PhaseImpl(const T& val) noexcept {
    multiply(val);
  }

  template <class T>
  void multiply(const T& val) noexcept {
    if (val == 0)
      makeNull();
    else if (val < 0)
      sign_ *= -1;
  }

  template <class T>
  void divide(const T& val) noexcept {
    if (val == 0)
      makeNull();
    else if (val < 0)
      sign_ *= -1;
  }

  void flip() noexcept {
    sign_ *= -1;
  }

  void reset() noexcept {
    sign_ = 1;
  }

  int8_t getSign() const noexcept {
    return sign_;
  }

  void makeNull() noexcept {
    sign_ = 0;
  }

  bool isNull() const noexcept {
    return sign_ == 0;
  };

  operator std::int8_t() const noexcept {
    return sign_;
  }

private:
  std::int8_t sign_ = 1;
};

template <class T>
using Phase = PhaseImpl<dca::util::IsComplex_t<T>::value>;

  template<typename T, typename = bool>
struct IsPhase : public std::false_type {};

template<typename T>
struct IsPhase<T, typename std::enable_if_t<(std::is_same<T, PhaseImpl<true>>::value || std::is_same<T, PhaseImpl<false>>::value), bool>>
 : public std::true_type {};    
  
}  // namespace dca::math

#endif  // DCA_MATH_UTIL_PHASE_HPP
