// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements the atomic convolution for NFFT.

#ifndef DCA_MATH_NFFT_NFFT_ATOMIC_CONVOLUTION_HPP
#define DCA_MATH_NFFT_NFFT_ATOMIC_CONVOLUTION_HPP

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

template <int max_count, int count>
struct NfftAtomicConvolutionImpl {
  template <typename ScalarType>
  inline static void execute_linear(ScalarType* f, const ScalarType* M, const ScalarType* y) {
    f[count] += (M[0 + 2 * count] * y[0] + M[1 + 2 * count] * y[1]);

    NfftAtomicConvolutionImpl<max_count, count + 1>::execute_linear(f, M, y);
  }

  template <typename ScalarType>
  inline static void execute_cubic(ScalarType* f, const ScalarType* M, const ScalarType* y) {
    f[count] += (M[0 + 4 * count] * y[0] + M[1 + 4 * count] * y[1] + M[2 + 4 * count] * y[2] +
                 M[3 + 4 * count] * y[3]);

    NfftAtomicConvolutionImpl<max_count, count + 1>::execute_cubic(f, M, y);
  }

  template <typename ScalarType>
  inline static void execute_M_y_2(ScalarType* f, const ScalarType* M, const int LDA,
                                   const ScalarType* y) {
    execute_Mt_y_1(f, M + 0 * LDA, y + 0);
    execute_Mt_y_1(f, M + 1 * LDA, y + 1);
  }

  template <typename ScalarType>
  inline static void execute_M_y_4(ScalarType* f, const ScalarType* M, const int LDA,
                                   const ScalarType* y) {
    execute_Mt_y_1(f, M + 0 * LDA, y + 0);
    execute_Mt_y_1(f, M + 1 * LDA, y + 1);
    execute_Mt_y_1(f, M + 0 * LDA, y + 2);
    execute_Mt_y_1(f, M + 1 * LDA, y + 3);
  }

  template <typename ScalarType>
  inline static void execute_Mt_y_1(ScalarType* f, const ScalarType* M, const ScalarType* y) {
    f[count] += M[count] * y[0];

    NfftAtomicConvolutionImpl<max_count, count + 1>::execute_Mt_y_1(f, M, y);
  }

  template <typename ScalarType>
  inline static void execute_Mt_y_2(ScalarType* f, const ScalarType* M, const ScalarType* y) {
    f[count] += (M[0 + 2 * count] * y[0] + M[1 + 2 * count] * y[1]);

    NfftAtomicConvolutionImpl<max_count, count + 1>::execute_Mt_y_2(f, M, y);
  }

  template <typename ScalarType>
  inline static void execute_Mt_y_4(ScalarType* f, const ScalarType* M, const ScalarType* y) {
    f[count] += (M[0 + 4 * count] * y[0] + M[1 + 4 * count] * y[1] + M[2 + 4 * count] * y[2] +
                 M[3 + 4 * count] * y[3]);

    NfftAtomicConvolutionImpl<max_count, count + 1>::execute_Mt_y_4(f, M, y);
  }
};

template <int max_count>
struct NfftAtomicConvolutionImpl<max_count, max_count> {
  template <typename ScalarType>
  inline static void execute_linear(ScalarType* /*f*/, const ScalarType* /*y*/,
                                    const ScalarType* /*M*/) {}

  template <typename ScalarType>
  inline static void execute_cubic(ScalarType* /*f*/, const ScalarType* /*y*/,
                                   const ScalarType* /*M*/) {}

  template <typename ScalarType>
  inline static void execute_Mt_y_2(ScalarType* /*f*/, const ScalarType* /*M*/,
                                    const ScalarType* /*y*/) {}

  template <typename ScalarType>
  inline static void execute_Mt_y_4(ScalarType* /*f*/, const ScalarType* /*M*/,
                                    const ScalarType* /*y*/) {}
};

template <int oversampling>
struct NfftAtomicConvolution {
  template <typename ScalarType>
  inline static void execute_linear(ScalarType* f, const ScalarType* y, const ScalarType* M) {
    static_assert(oversampling > 1, "Invalid oversampling size.");
    NfftAtomicConvolutionImpl<2 * oversampling + 1, 1>::execute_linear(f, y, M);
  }

  template <typename ScalarType>
  inline static void execute_cubic(ScalarType* f, const ScalarType* y, const ScalarType* M) {
    static_assert(oversampling > 1, "Invalid oversampling size.");
    NfftAtomicConvolutionImpl<2 * oversampling + 1, 1>::execute_cubic(f, y, M);
  }
};

}  // namespace nfft
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_NFFT_NFFT_ATOMIC_CONVOLUTION_HPP
