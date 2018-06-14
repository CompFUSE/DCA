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

template <typename scalartype>
class akima_interpolation {
public:
  akima_interpolation(int n);
  ~akima_interpolation();

  void initialize(scalartype* x, scalartype* y);
  void initialize_periodic(scalartype* x, scalartype* y);

  scalartype evaluate(scalartype x);

  scalartype get_alpha(int l, int i);

private:
  void compute_coefficients(scalartype* x_array, scalartype* y_array);

  int size;

  scalartype* X;
  scalartype* Y;

  scalartype* a;
  scalartype* b;
  scalartype* c;
  scalartype* d;

  scalartype* m;
  scalartype* _m;
};

template <typename scalartype>
akima_interpolation<scalartype>::akima_interpolation(int n)
    : size(n),

      X(nullptr),
      Y(nullptr),

      a(nullptr),
      b(nullptr),
      c(nullptr),
      d(nullptr),

      m(nullptr),
      _m(nullptr) {
  X = new scalartype[size];
  Y = new scalartype[size];

  a = new scalartype[size];
  b = new scalartype[size];
  c = new scalartype[size];
  d = new scalartype[size];

  _m = new scalartype[size + 4];
}

template <typename scalartype>
akima_interpolation<scalartype>::~akima_interpolation() {
  if (X != nullptr)
    delete[] X;
  if (Y != nullptr)
    delete[] Y;

  if (a != nullptr)
    delete[] a;
  if (b != nullptr)
    delete[] b;
  if (c != nullptr)
    delete[] c;
  if (d != nullptr)
    delete[] d;

  if (_m != nullptr)
    delete[] _m;
}

template <typename scalartype>
void akima_interpolation<scalartype>::initialize(scalartype* x_array, scalartype* y_array) {
  for (int i = 0; i < size; i++) {
    X[i] = x_array[i];
    Y[i] = y_array[i];
  }

  m = _m + 2; /* offset so we can address the -1,-2 components */

  for (int i = 0; i <= size - 2; i++)
    m[i] = (y_array[i + 1] - y_array[i]) / (x_array[i + 1] - x_array[i]);

  /* non-periodic boundary conditions */
  m[-2] = 3.0 * m[0] - 2.0 * m[1];
  m[-1] = 2.0 * m[0] - m[1];
  m[size - 1] = 2.0 * m[size - 2] - m[size - 3];
  m[size] = 3.0 * m[size - 2] - 2.0 * m[size - 3];

  compute_coefficients(x_array, y_array);
}

template <typename scalartype>
void akima_interpolation<scalartype>::initialize_periodic(scalartype* x_array, scalartype* y_array) {
  m = _m + 2; /* offset so we can address the -1,-2 components */

  for (int i = 0; i <= size - 2; i++)
    m[i] = (y_array[i + 1] - y_array[i]) / (x_array[i + 1] - x_array[i]);

  /* periodic boundary conditions */
  m[-2] = m[size - 1 - 2];
  m[-1] = m[size - 1 - 1];
  m[size - 1] = m[0];
  m[size] = m[1];

  compute_coefficients(x_array, y_array);
}

template <typename scalartype>
scalartype akima_interpolation<scalartype>::evaluate(scalartype x) {
  if (x < X[0] + 1.e-6)
    return Y[0];

  if (x >= X[size - 1] - 1.e-6)
    return Y[size - 1];

  int index = -1;
  for (int i = 0; i < size - 1; ++i)
    if (X[i] <= x and x < X[i + 1])
      index = i;

  assert(index > -1 and index < size);

  const scalartype x_lo = X[index];
  const scalartype delx = x - x_lo;

  const scalartype a0 = a[index];
  const scalartype a1 = b[index];
  const scalartype a2 = c[index];
  const scalartype a3 = d[index];

  return a0 + delx * (a1 + delx * (a2 + a3 * delx));
}

template <typename scalartype>
scalartype akima_interpolation<scalartype>::get_alpha(int l, int i) {
  assert(i > -1 and i < size - 1);

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

template <typename scalartype>
void akima_interpolation<scalartype>::compute_coefficients(scalartype* x_array, scalartype* y_array) {
  for (int i = 0; i < (size - 1); i++) {
    a[i] = y_array[i];

    const scalartype NE = std::abs(m[i + 1] - m[i]) + fabs(m[i - 1] - m[i - 2]);

    if (NE == 0.0) {
      b[i] = m[i];
      c[i] = 0.0;
      d[i] = 0.0;
    }
    else {
      const scalartype h_i = x_array[i + 1] - x_array[i];
      const scalartype NE_next = std::abs(m[i + 2] - m[i + 1]) + fabs(m[i] - m[i - 1]);
      const scalartype alpha_i = std::abs(m[i - 1] - m[i - 2]) / NE;

      scalartype alpha_ip1;
      scalartype tL_ip1;

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

}  // interpolation
}  // math
}  // dca

#endif  // DCA_MATH_INTERPOLATION_AKIMA_INTERPOLATION_HPP
