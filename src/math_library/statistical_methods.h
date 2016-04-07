//-*-C++-*-

#ifndef MATH_LIBRARY_STATISTICAL_METHODS_H
#define MATH_LIBRARY_STATISTICAL_METHODS_H

namespace math_algorithms {

template <typename scalartype>
class statistical_methods {
public:
  static scalartype mean(std::vector<scalartype>& vec);
  static scalartype variance(std::vector<scalartype>& vec);
  static scalartype standard_deviation(std::vector<scalartype>& vec);

  static std::vector<scalartype> auto_correlation_time(std::vector<scalartype>& vec);
  static std::vector<scalartype> auto_correlation_time_2(std::vector<scalartype>& vec);
};

template <typename scalartype>
scalartype statistical_methods<scalartype>::mean(std::vector<scalartype>& vec) {
  assert(int(vec.size()) > 0);

  scalartype result(0);

  for (size_t i = 0; i < vec.size(); i++)
    result += vec[i];

  return result / scalartype(vec.size());
}

template <typename scalartype>
scalartype statistical_methods<scalartype>::variance(std::vector<scalartype>& vec) {
  assert(int(vec.size()) > 1);

  scalartype m = mean(vec);
  scalartype result(0);

  for (size_t i = 0; i < vec.size(); i++)
    result += square((vec[i] - m));

  return result / (scalartype(vec.size() - 1));
}

template <typename scalartype>
scalartype statistical_methods<scalartype>::standard_deviation(std::vector<scalartype>& vec) {
  return sqrt(variance(vec));
}

template <typename scalartype>
std::vector<scalartype> statistical_methods<scalartype>::auto_correlation_time(
    std::vector<scalartype>& vec) {
  int N = vec.size();

  std::vector<scalartype> result(0);

  for (int l = 1; l < 2000; l++) {
    int bin_size = l;  // pow(2.,l);

    std::vector<scalartype> mean_vec;

    for (int z = 0; z < N / bin_size; z++) {
      if (bin_size - 1 + z * bin_size < N) {
        scalartype mean = 0;

        for (int i = 0; i < bin_size; i++)
          mean += vec[i + z * bin_size] / double(bin_size);

        mean_vec.push_back(mean);
      }
    }

    if (int(mean_vec.size()) > 1)
      result.push_back(variance(mean_vec));
  }

  scalartype num = result[0];
  for (size_t l = 0; l < result.size(); l++)
    result[l] /= num;

  return result;
}

template <typename scalartype>
std::vector<scalartype> statistical_methods<scalartype>::auto_correlation_time_2(
    std::vector<scalartype>& vec) {
  int N = vec.size();

  scalartype m = 0;
  for (int i = 0; i < N; i++)
    m += vec[i] / double(N);

  std::vector<scalartype> result(N / 2, 0.);

  for (int t = 0; t < N / 2; t++)
    for (int i = 0; i < N - t; i++)
      result[t] += (vec[i] - m) * (vec[i + t] - m) / double(N - t);

  return result;
}
}
#endif  // MATH_LIBRARY_STATISTICAL_METHODS_H
