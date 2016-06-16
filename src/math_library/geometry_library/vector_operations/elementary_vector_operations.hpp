// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description: TBA
//
// TODO: - Remove duplicated code.
//       - Check function bodies for const correctness.

#ifndef MATH_LIBRARY_GEOMETRY_LIBRARY_VECTOR_OPERATIONS_ELEMENTARY_VECTOR_OPERATIONS_HPP
#define MATH_LIBRARY_GEOMETRY_LIBRARY_VECTOR_OPERATIONS_ELEMENTARY_VECTOR_OPERATIONS_HPP

#include <algorithm>
#include "comp_library/linalg/linalg.hpp"

namespace VECTOR_OPERATIONS {
template <typename scalartype>
inline scalartype NORM(const std::vector<scalartype>& v1) {
  scalartype norm = 0;
  for (size_t i = 0; i < v1.size(); ++i) {
    norm += v1[i] * v1[i];
  }

  return std::sqrt(norm);
}

template <typename scalartype>
inline scalartype NORM(const std::vector<scalartype>& v1, const std::vector<scalartype>& v2) {
  assert(v1.size() == v2.size());

  scalartype norm = 0;
  for (size_t i = 0; i < v1.size(); ++i) {
    norm += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  }

  return std::sqrt(norm);
}

template <typename scalartype>
inline scalartype L2_NORM(const std::vector<scalartype>& v1) {
  scalartype norm = 0;
  for (size_t i = 0; i < v1.size(); ++i) {
    norm += v1[i] * v1[i];
  }

  return norm;
}

template <typename scalartype>
inline scalartype L2_NORM(const std::vector<scalartype>& v1, const std::vector<scalartype>& v2) {
  assert(v1.size() == v2.size());

  scalartype norm = 0;
  for (size_t i = 0; i < v1.size(); ++i) {
    norm += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  }

  return norm;
}

template <typename scalartype>
inline scalartype L2_NORM(const int N, const scalartype* v1, const scalartype* v2) {
  scalartype norm = 0;
  for (int i = 0; i < N; ++i) {
    norm += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  }

  return norm;
}

bool SAME_VECTOR(const std::vector<double> v1, const std::vector<double> v2) {
  return (L2_NORM(v1, v2) < 1.e-6);
}

template <typename scalartype>
bool IS_LARGER_VECTOR(const std::vector<scalartype>& v1, const std::vector<scalartype>& v2) {
  assert(v1.size() == v2.size());

  for (size_t i = 0; i < v1.size(); ++i) {
    if (std::fabs(v1[i] - v2[i]) > 1.e-6) {
      return v1[i] < v2[i];
    }
  }

  return false;
}

template <typename scalartype>
bool HAS_LARGER_NORM(const std::vector<scalartype> v1, const std::vector<scalartype> v2) {
  assert(v1.size() == v2.size());

  if (std::abs(L2_NORM(v1) - L2_NORM(v2)) < 1.e-6) {
    return IS_LARGER_VECTOR(v1, v2);
  }

  if (L2_NORM(v1) < L2_NORM(v2)) {
    return true;
  }
  else {
    return false;
  }
}

template <typename scalartype>
inline std::vector<scalartype> SCALE(const scalartype a, const std::vector<scalartype>& v) {
  std::vector<scalartype> tmp = v;

  for (size_t i = 0; i < tmp.size(); ++i) {
    tmp[i] *= a;
  }

  return tmp;
}

template <typename scalartype>
inline std::vector<scalartype> ADD(const std::vector<scalartype> v1,
                                   const std::vector<scalartype> v2) {
  assert(v1.size() == v2.size());

  std::vector<scalartype> tmp(v1.size());

  for (size_t i = 0; i < v1.size(); ++i) {
    tmp[i] = v1[i] + v2[i];
  }

  return tmp;
}

template <typename scalartype>
inline std::vector<scalartype> SUBTRACT(const std::vector<scalartype> v1,
                                        const std::vector<scalartype> v2) {
  assert(v1.size() == v2.size());

  std::vector<scalartype> tmp(v1.size());

  for (size_t i = 0; i < v1.size(); ++i) {
    tmp[i] = v2[i] - v1[i];
  }

  return tmp;
}

template <typename scalartype>
inline scalartype DOT_PRODUCT(const std::vector<scalartype>& v1, const std::vector<scalartype>& v2) {
  assert(v1.size() == v2.size());

  scalartype norm = 0;
  for (size_t i = 0; i < v1.size(); ++i) {
    norm += v1[i] * v2[i];
  }

  return norm;
}

template <typename scalartype>
inline void CROSS_PRODUCT(const std::vector<scalartype>& v1, const std::vector<scalartype>& v2,
                          std::vector<scalartype>& v3) {
  assert(v1.size() == 3);
  assert(v2.size() == 3);
  assert(v3.size() == 3);

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = -(v1[0] * v2[2] - v1[2] * v2[0]);
  v3[1] = v1[0] * v2[1] - v1[1] * v2[0];
}

template <typename scalartype>
inline std::vector<scalartype> COORDINATES(const std::vector<scalartype>& r,
                                           const std::vector<std::vector<scalartype>>& B) {
  assert(r.size() == B.size());

  int N = r.size();

  std::vector<scalartype> basis(N * N);  // scalartype basis[N*N];

  for (int d1 = 0; d1 < N; ++d1) {
    for (int d0 = 0; d0 < N; ++d0) {
      basis[d0 + d1 * N] = B[d1][d0];
    }
  }

  std::vector<scalartype> coordinate = r;
  LIN_ALG::GESV<LIN_ALG::CPU>::execute(N, &basis[0], &coordinate[0]);

  return coordinate;
}

template <typename scalartype>
inline void COORDINATES(const std::vector<scalartype>& v1, const std::vector<scalartype>& v2,
                        const std::vector<scalartype>& vec, std::vector<scalartype>& coor) {
  assert(v1.size() == 2);
  assert(v2.size() == 2);
  assert(vec.size() == 2);
  assert(coor.size() == 2);

  int N = 2;

  std::vector<scalartype> basis(N * N);  // scalartype basis[N*N];

  basis[0 + 2 * 0] = v1[0];
  basis[1 + 2 * 0] = v1[1];

  basis[0 + 2 * 1] = v2[0];
  basis[1 + 2 * 1] = v2[1];

  LIN_ALG::GESV<LIN_ALG::CPU>::execute(N, &basis[0], &coor[0]);
}

template <typename scalartype>
inline void COORDINATES(const std::vector<scalartype>& v1, const std::vector<scalartype>& v2,
                        const std::vector<scalartype>& v3, const std::vector<scalartype>& vec,
                        std::vector<scalartype>& coor) {
  assert(v1.size() == 3);
  assert(v2.size() == 3);
  assert(v3.size() == 3);
  assert(vec.size() == 3);
  assert(coor.size() == 3);

  int N = 3;

  std::vector<scalartype> basis(N * N);  // scalartype basis[N*N];

  basis[0 + 3 * 0] = v1[0];
  basis[1 + 3 * 0] = v1[1];
  basis[2 + 3 * 0] = v1[2];

  basis[0 + 3 * 1] = v2[0];
  basis[1 + 3 * 1] = v2[1];
  basis[2 + 3 * 1] = v2[2];

  basis[0 + 3 * 2] = v3[0];
  basis[1 + 3 * 2] = v3[1];
  basis[2 + 3 * 2] = v3[2];

  LIN_ALG::GESV<LIN_ALG::CPU>::execute(N, &basis[0], &coor[0]);
}

template <typename scalartype>
inline scalartype VOLUME(const std::vector<scalartype>& v1, const std::vector<scalartype>& v2) {
  assert(v1.size() == 2);
  assert(v2.size() == 2);

  return std::abs(v1[0] * v2[1] - v1[1] * v2[0]);
}

template <typename scalartype>
inline scalartype VOLUME(const std::vector<scalartype>& v1, const std::vector<scalartype>& v2,
                         const std::vector<scalartype>& v3) {
  assert(v1.size() == 3);
  assert(v2.size() == 3);
  assert(v3.size() == 3);

  return std::fabs(v1[0] * v2[1] * v3[2] + v2[0] * v3[1] * v1[2] + v1[1] * v2[2] * v3[0] -
                   v1[2] * v2[1] * v3[0] - v2[0] * v1[1] * v3[2] - v2[2] * v3[1] * v1[0]);
}

template <typename scalartype>
inline void PRINT(const std::vector<scalartype> v1) {
  std::cout << std::scientific;
  std::cout.precision(6);

  for (size_t i = 0; i < v1.size(); ++i) {
    std::cout << v1[i] << "\t";
  }
}

template <typename scalartype1, typename scalartype2>
inline scalartype1 MINIMIZE(const int N, const scalartype1* x, const scalartype2* y) {
  int index = std::min_element(y, y + N) - y;
  assert(index > -1 && index < N);

  if (index == 0 || index == N - 1) {
    return x[index];
  }
  else {
    scalartype1* x_ptr = x + index - 1;
    scalartype1* y_ptr = y + index - 1;

    scalartype1 a = -((-(x_ptr[1] * y_ptr[0]) + x_ptr[2] * y_ptr[0] + x_ptr[0] * y_ptr[1] -
                       x_ptr[2] * y_ptr[1] - x_ptr[0] * y_ptr[2] + x_ptr[1] * y_ptr[2]) /
                      ((x_ptr[1] - x_ptr[2]) * (std::pow(x_ptr[0], 2) - x_ptr[0] * x_ptr[1] -
                                                x_ptr[0] * x_ptr[2] + x_ptr[1] * x_ptr[2])));
    scalartype1 b = -((std::pow(x_ptr[1], 2) * y_ptr[0] - std::pow(x_ptr[2], 2) * y_ptr[0] -
                       std::pow(x_ptr[0], 2) * y_ptr[1] + std::pow(x_ptr[2], 2) * y_ptr[1] +
                       std::pow(x_ptr[0], 2) * y_ptr[2] - std::pow(x_ptr[1], 2) * y_ptr[2]) /
                      ((x_ptr[0] - x_ptr[1]) * (x_ptr[0] - x_ptr[2]) * (x_ptr[1] - x_ptr[2])));

    return -b / (2. * a);
  }
}

template <typename scalartype>
inline scalartype MINIMIZE(const std::vector<scalartype>& x, const std::vector<scalartype>& y) {
  assert(x.size() == y.size());
  return MINIMIZE(int(x.size()), &x[0], &y[0]);
}

template <typename scalartype1, typename scalartype2>
inline scalartype1 MAXIMIZE(const int N, const scalartype1* x, const scalartype2* y) {
  int index = std::max_element(y, y + N) - y;
  assert(index > -1 && index < N);

  if (index == 0 || index == N - 1) {
    return x[index];
  }
  else {
    scalartype1* x_ptr = x + index - 1;
    scalartype1* y_ptr = y + index - 1;

    scalartype1 a = -((-(x_ptr[1] * y_ptr[0]) + x_ptr[2] * y_ptr[0] + x_ptr[0] * y_ptr[1] -
                       x_ptr[2] * y_ptr[1] - x_ptr[0] * y_ptr[2] + x_ptr[1] * y_ptr[2]) /
                      ((x_ptr[1] - x_ptr[2]) * (std::pow(x_ptr[0], 2) - x_ptr[0] * x_ptr[1] -
                                                x_ptr[0] * x_ptr[2] + x_ptr[1] * x_ptr[2])));
    scalartype1 b = -((std::pow(x_ptr[1], 2) * y_ptr[0] - std::pow(x_ptr[2], 2) * y_ptr[0] -
                       std::pow(x_ptr[0], 2) * y_ptr[1] + std::pow(x_ptr[2], 2) * y_ptr[1] +
                       std::pow(x_ptr[0], 2) * y_ptr[2] - std::pow(x_ptr[1], 2) * y_ptr[2]) /
                      ((x_ptr[0] - x_ptr[1]) * (x_ptr[0] - x_ptr[2]) * (x_ptr[1] - x_ptr[2])));

    return -b / (2. * a);
  }
}
}  // VECTOR_OPERATIONS

#endif  // MATH_LIBRARY_GEOMETRY_LIBRARY_VECTOR_OPERATIONS_ELEMENTARY_VECTOR_OPERATIONS_HPP
