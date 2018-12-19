// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides utility functions to do various vector operations.

#ifndef DCA_MATH_UTIL_VECTOR_OPERATIONS_HPP
#define DCA_MATH_UTIL_VECTOR_OPERATIONS_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <type_traits>
#include <vector>

#include "dca/linalg/lapack/solve.hpp"

namespace dca {
namespace math {
namespace util {
// dca::math::util::

// Prints the elements of the vector v.
template <typename T>
void print(const std::vector<T>& v) {
  std::cout << std::scientific;
  std::cout.precision(6);

  for (std::size_t i = 0; i < v.size(); ++i)
    std::cout << v[i] << "\t";
}

// Scales the vector v by the scalar a, i.e. computes and returns a new vector w with elements
// w[i] = a * v[i], 0 <= i < v.size().
template <typename T>
std::vector<T> scale(const T a, const std::vector<T>& v) {
  std::vector<T> w = v;

  for (std::size_t i = 0; i < w.size(); ++i)
    w[i] *= a;

  return w;
}

// Computes and returns the sum of two vectors x and y.
// Preconditions: x.size() == y.size().
template <typename T>
std::vector<T> add(const std::vector<T>& x, const std::vector<T>& y) {
  assert(x.size() == y.size());

  std::vector<T> res(x);

  for (std::size_t i = 0; i < res.size(); ++i)
    res[i] += y[i];

  return res;
}

// Subtracts the vector x from the vector y and returns the result, i.e. the vector z with elements
// z[i] = y[i] - x[i], 0 <= i < x.size().
// Preconditions: x.size() == y.size().
template <typename T>
std::vector<T> subtract(const std::vector<T>& x, const std::vector<T>& y) {
  assert(x.size() == y.size());

  std::vector<T> res(y);

  for (std::size_t i = 0; i < res.size(); ++i)
    res[i] -= x[i];

  return res;
}

// Computes the inner product of two vectors x and y.
// Preconditions: x.size() == y.size().
template <typename T>
T innerProduct(const std::vector<T>& x, const std::vector<T>& y) {
  assert(x.size() == y.size());

  T res = 0;
  for (std::size_t i = 0; i < x.size(); ++i)
    res += x[i] * y[i];

  return res;
}
template <typename T>
std::complex<T> innerProduct(const std::vector<std::complex<T>>& x,
                             const std::vector<std::complex<T>>& y) {
  assert(x.size() == y.size());

  std::complex<T> res(0.);
  for (std::size_t i = 0; i < x.size(); ++i)
    res += x[i] * std::conj(y[i]);

  return res;
}

// Treats scalars as vectors of size 1.
template <typename T>
T innerProduct(const T x, const T y) {
  return x * y;
}
template <typename T>
std::complex<T> innerProduct(const std::complex<T> x, const std::complex<T> y) {
  return x * std::conj(y);
}

// Computes the square of the L^2 norm of the vector x.
template <typename T>
auto l2Norm2(const std::vector<T>& x) {
  return std::real(innerProduct(x, x));
}

// Computes the L^2 norm of the vector x.
template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, T> l2Norm(const std::vector<T>& x) {
  return std::sqrt(l2Norm2(x));
}
template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, T> l2Norm(const std::vector<std::complex<T>>& x) {
  return std::sqrt(l2Norm2(x));
}

// Computes the square of the distance of two vectors x and y in terms of the L^2 norm of their
// difference.
// Preconditions: x.size() == y.size().
template <typename T>
auto distance2(const std::vector<T>& x, const std::vector<T>& y) {
  return l2Norm2(subtract(y, x));
}

// Computes the distance of two vectors x and y in terms of the L^2 norm of their
// difference.
// Preconditions: x.size() == y.size().
template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, T> distance(const std::vector<T>& x,
                                                               const std::vector<T>& y) {
  return std::sqrt(distance2(x, y));
}
template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, T> distance(const std::vector<std::complex<T>>& x,
                                                               const std::vector<std::complex<T>>& y) {
  return std::sqrt(distance2(x, y));
}

// Computes the area of the parallelogram that the two 2D vectors x and y span.
// Preconditions: x.size() == y.size() == 2.
template <typename T>
std::enable_if_t<std::is_arithmetic<T>::value, T> area(const std::vector<T>& x,
                                                       const std::vector<T>& y) {
  assert(x.size() == 2);
  assert(y.size() == 2);
  return std::abs(x[0] * y[1] - x[1] * y[0]);
}

// Computes the volume of the parallelpiped that the three 3D vectors x, y and z span.
// Preconditions: x.size() == y.size() == z.size() == 3.
template <typename T>
std::enable_if_t<std::is_arithmetic<T>::value, T> volume(const std::vector<T>& x,
                                                         const std::vector<T>& y,
                                                         const std::vector<T>& z) {
  assert(x.size() == 3);
  assert(y.size() == 3);
  assert(z.size() == 3);

  return std::abs(x[0] * y[1] * z[2] + y[0] * z[1] * x[2] + x[1] * y[2] * z[0] -
                  x[2] * y[1] * z[0] - y[0] * x[1] * z[2] - y[2] * z[1] * x[0]);
}

// Computes the coordinates of the vector r in the basis defined by the elements of b.
// Preconditions: - r.size() == b.size()
//                - r.size() == b[i].size(), 0 <= i < r.size().
template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, std::vector<T>> coordinates(
    const std::vector<T>& r, const std::vector<std::vector<T>>& b) {
  const std::size_t n = r.size();

  assert(b.size() == n);
  for (std::size_t i = 0; i < b.size(); ++i)
    assert(b[i].size() == n);

  std::vector<T> basis(n * n);

  for (std::size_t d1 = 0; d1 < r.size(); ++d1) {
    for (std::size_t d0 = 0; d0 < r.size(); ++d0) {
      basis[d0 + d1 * n] = b[d1][d0];
    }
  }

  std::vector<T> coord = r;
  linalg::lapack::solve(n, &basis[0], n, &coord[0]);

  return coord;
}

// Compare function for two vectors x and y.
// Returns the result of x[i] < y[i] for the first index i with |x[i] - y[i]| > 1.e-6.
// Preconditions: x.size() == y.size().
template <typename T>
std::enable_if_t<std::is_arithmetic<T>::value, bool> isLessVector(const std::vector<T>& x,
                                                                  const std::vector<T>& y) {
  const double tol = 1.e-6;
  assert(x.size() == y.size());

  for (std::size_t i = 0; i < x.size(); ++i) {
    if (std::abs(x[i] - y[i]) > tol) {
      return x[i] < y[i];
    }
  }

  return false;
}

// Compare function for two vectors x and y.
// Returns true, if the L^2 norm of x is smaller than the L^2 norm of y.
// If the squares of the norms differ by less than 1.e-6, the function falls back on isLessVector.
// Preconditions: x.size() == y.size().
template <typename T>
std::enable_if_t<std::is_arithmetic<T>::value, bool> hasSmallerNorm(const std::vector<T>& x,
                                                                    const std::vector<T>& y) {
  assert(x.size() == y.size());

  const double tol = 1.e-6;

  if (std::abs(l2Norm2(x) - l2Norm2(y)) < tol)
    return isLessVector(x, y);

  else
    return l2Norm2(x) < l2Norm(y);
}

// Compare function for two vectors x and y.
// Returns true, if the squared distance between the vectors is smaller than 1.e-6.
template <typename T>
std::enable_if_t<std::is_arithmetic<T>::value, bool> isSameVector(const std::vector<T>& x,
                                                                  const std::vector<T>& y) {
  const double tol = 1.e-6;
  return (distance2(x, y) < tol);
}

// Compare two vectors element-wise.
template <typename T>
bool operator==(const std::vector<T>& x, const std::vector<T>& y) {
  if (x.size() != y.size())
    return false;
  for (int i = 0; i < x.size(); ++i)
    if (!(x[i] == y[i]))
      return false;
  return true;
}

}  // util
}  // math
}  // dca

#endif  // DCA_MATH_UTIL_VECTOR_OPERATIONS_HPP
