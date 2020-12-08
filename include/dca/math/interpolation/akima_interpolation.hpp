// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements Akima's interpolation method.

#ifndef DCA_MATH_INTERPOLATION_AKIMA_INTERPOLATION_HPP
#define DCA_MATH_INTERPOLATION_AKIMA_INTERPOLATION_HPP

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace dca {
namespace math {
namespace interpolation {
// dca::math::interpolation::

template <typename Scalar>
class akima_interpolation {
public:
  akima_interpolation(size_t n);
  ~akima_interpolation() = default;

  void initialize(const std::vector<Scalar>& x, const std::vector<Scalar>& y);
  void initialize_periodic(const std::vector<Scalar>& x, const std::vector<Scalar>& y);

  Scalar evaluate(const Scalar x) const;

  Scalar get_alpha(int l, int i) const;

  std::size_t size() const {
    return size_;
  }

private:
  void compute_coefficients(const std::vector<Scalar>& x, const std::vector<Scalar>& y);

  std::size_t size_;

  std::vector<Scalar> X;
  std::vector<Scalar> Y;

  std::vector<Scalar> a;
  std::vector<Scalar> b;
  std::vector<Scalar> c;
  std::vector<Scalar> d;

  std::vector<Scalar> m_data;
  Scalar* m;
};

template <typename Scalar>
akima_interpolation<Scalar>::akima_interpolation(size_t n)
    : size_(n),

      X(n),
      Y(n),

      a(n),
      b(n),
      c(n),
      d(n),

      m_data(n + 4) {
  m = m_data.data() + 2; /* offset so we can address the -1,-2 components */
}

template <typename Scalar>
void akima_interpolation<Scalar>::initialize(const std::vector<Scalar>& x_array,
                                             const std::vector<Scalar>& y_array) {
  assert(size() == x_array.size());
  assert(size() == y_array.size());

  for (int i = 0; i < size_; i++) {
    X[i] = x_array[i];
    Y[i] = y_array[i];
  }

  for (int i = 0; i <= size_ - 2; i++)
    m[i] = (y_array[i + 1] - y_array[i]) / (x_array[i + 1] - x_array[i]);

  /* non-periodic boundary conditions */
  m[-2] = 3.0 * m[0] - 2.0 * m[1];
  m[-1] = 2.0 * m[0] - m[1];
  m[size_ - 1] = 2.0 * m[size_ - 2] - m[size_ - 3];
  m[size_] = 3.0 * m[size_ - 2] - 2.0 * m[size_ - 3];

  compute_coefficients(x_array, y_array);
}

template <typename Scalar>
void akima_interpolation<Scalar>::initialize_periodic(const std::vector<Scalar>& x_array, const std::vector<Scalar>& y_array) {
  assert(size() == x_array.size());
  assert(size() == y_array.size());

  for (int i = 0; i <= size_ - 2; i++)
    m[i] = (y_array[i + 1] - y_array[i]) / (x_array[i + 1] - x_array[i]);

  /* periodic boundary conditions */
  m[-2] = m[size_ - 1 - 2];
  m[-1] = m[size_ - 1 - 1];
  m[size_ - 1] = m[0];
  m[size_] = m[1];

  compute_coefficients(x_array, y_array);
}

template <typename Scalar>
Scalar akima_interpolation<Scalar>::evaluate(Scalar x) const {
  if (x < X[0] + 1.e-6)
    return Y[0];

  if (x >= X[size_ - 1] - 1.e-6)
    return Y[size_ - 1];

  int index = -1;
  for (int i = 0; i < size_ - 1; ++i)
    if (X[i] <= x and x < X[i + 1])
      index = i;

  assert(index > -1 and index < size_);

  const Scalar x_lo = X[index];
  const Scalar delx = x - x_lo;

  const Scalar a0 = a[index];
  const Scalar a1 = b[index];
  const Scalar a2 = c[index];
  const Scalar a3 = d[index];

  return a0 + delx * (a1 + delx * (a2 + a3 * delx));
}

template <typename Scalar>
Scalar akima_interpolation<Scalar>::get_alpha(int l, int i) const {
  assert(i > -1 and i < size_ - 1);

  switch (l) {
    case 0:
      return a[i];

    case 1:
      return b[i];

    case 2:
      return c[i];

    case 3:
      return d[i];

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename Scalar>
void akima_interpolation<Scalar>::compute_coefficients(const std::vector<Scalar>& x_array, const std::vector<Scalar>& y_array) {
  for (int i = 0; i < (size_ - 1); i++) {
    a[i] = y_array[i];

    const Scalar NE = std::abs(m[i + 1] - m[i]) + fabs(m[i - 1] - m[i - 2]);

    if (NE == 0.0) {
      b[i] = m[i];
      c[i] = 0.0;
      d[i] = 0.0;
    }
    else {
      const Scalar h_i = x_array[i + 1] - x_array[i];
      const Scalar NE_next = std::abs(m[i + 2] - m[i + 1]) + fabs(m[i] - m[i - 1]);
      const Scalar alpha_i = std::abs(m[i - 1] - m[i - 2]) / NE;

      Scalar alpha_ip1;
      Scalar tL_ip1;

      if (NE_next == 0.0) {
        tL_ip1 = m[i];
      }
      else {
        alpha_ip1 = std::abs(m[i] - m[i - 1]) / NE_next;
        tL_ip1 = (1.0 - alpha_ip1) * m[i] + alpha_ip1 * m[i + 1];
      }

      b[i] = (1.0 - alpha_i) * m[i - 1] + alpha_i * m[i];
      c[i] = (3.0 * m[i] - 2.0 * b[i] - tL_ip1) / h_i;
      d[i] = (b[i] + tL_ip1 - 2.0 * m[i]) / (h_i * h_i);
    }
  }
}

template <typename Scalar>
class akima_interpolation<std::complex<Scalar>> {
public:
  akima_interpolation(size_t n);
  ~akima_interpolation() = default;

  void initialize(const std::vector<Scalar>& x, const std::vector<std::complex<Scalar>>& y);

  std::complex<Scalar> evaluate(const Scalar x) const;

  std::complex<Scalar> get_alpha(int l, int i) const;

  std::size_t size() const {
    return interpolations_[0].size();
  }

private:
  std::array<akima_interpolation<Scalar>, 2> interpolations_;
};

template <typename Scalar>
akima_interpolation<std::complex<Scalar>>::akima_interpolation(size_t n) : interpolations_{n, n} {}

template <typename Scalar>
void akima_interpolation<std::complex<Scalar>>::initialize(const std::vector<Scalar>& x_array,
                                                           const std::vector<std::complex<Scalar>>& y_array) {
  std::array<std::vector<Scalar>, 2> y;
  const auto n = size();
  for (int re_im = 0; re_im < 2; ++re_im) {
    y[re_im].resize(n);
  }

  for (std::size_t i = 0; i < n; ++i) {
    y[0][i] = std::real(y_array[i]);
    y[1][i] = std::imag(y_array[i]);
  }

  for (int re_im = 0; re_im < 2; ++re_im)
    interpolations_[re_im].initialize(x_array, y[re_im]);
}

template <typename Scalar>
std::complex<Scalar> akima_interpolation<std::complex<Scalar>>::evaluate(Scalar x) const {
  return {interpolations_[0].evaluate(x), interpolations_[1].evaluate(x)};
}

template <typename Scalar>
std::complex<Scalar> akima_interpolation<std::complex<Scalar>>::get_alpha(int l, int i) const {
  return {interpolations_[0].get_alpha(l, i), interpolations_[1].get_alpha(l, i)};
}

}  // namespace interpolation
}  // namespace math
}  // namespace dca

#endif  // DCA_MATH_INTERPOLATION_AKIMA_INTERPOLATION_HPP
