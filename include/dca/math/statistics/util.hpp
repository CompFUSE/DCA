// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides utility methods to compute sample statistics.

#ifndef DCA_MATH_STATISTICS_UTIL_HPP
#define DCA_MATH_STATISTICS_UTIL_HPP

#include <cassert>
#include <cstdlib>
#include <type_traits>
#include <vector>

namespace dca {
namespace math {
namespace statistics {
namespace util {
// dca::math::statistics::util::

template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, T> mean(const std::vector<T>& vec) {
  assert(vec.size() > 0);

  T result(0);

  for (std::size_t i = 0; i < vec.size(); i++)
    result += vec[i];

  return result / T(vec.size());
}

template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, T> variance(const std::vector<T>& vec) {
  assert(vec.size() > 1);

  T m = mean(vec);
  T result(0);

  for (std::size_t i = 0; i < vec.size(); i++)
    result += (vec[i] - m) * (vec[i] - m);

  return result / (T(vec.size() - 1));
}

template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, T> standard_deviation(const std::vector<T>& vec) {
  return std::sqrt(variance(vec));
}

template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, std::vector<T>> auto_correlation_time(
    const std::vector<T>& vec) {
  std::size_t N = vec.size();

  std::vector<T> result(0);

  for (int l = 1; l < 2000; l++) {
    int bin_size = l;  // pow(2.,l);

    std::vector<T> mean_vec;

    for (int z = 0; z < N / bin_size; z++) {
      if (bin_size - 1 + z * bin_size < N) {
        T mean = 0;

        for (int i = 0; i < bin_size; i++)
          mean += vec[i + z * bin_size] / double(bin_size);

        mean_vec.push_back(mean);
      }
    }

    if (mean_vec.size() > 1)
      result.push_back(variance(mean_vec));
  }

  T num = result[0];
  for (std::size_t l = 0; l < result.size(); l++)
    result[l] /= num;

  return result;
}

template <typename T>
std::enable_if_t<std::is_floating_point<T>::value, std::vector<T>> auto_correlation_time_2(
    const std::vector<T>& vec) {
  std::size_t N = vec.size();

  T m = 0;
  for (int i = 0; i < N; i++)
    m += vec[i] / double(N);

  std::vector<T> result(N / 2, 0.);

  for (int t = 0; t < N / 2; t++)
    for (int i = 0; i < N - t; i++)
      result[t] += (vec[i] - m) * (vec[i + t] - m) / double(N - t);

  return result;
}

}  // util
}  // statistics
}  // math
}  // dca

#endif  // DCA_MATH_STATISTICS_UTIL_HPP
