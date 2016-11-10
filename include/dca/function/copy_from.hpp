// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements an efficient copy mechanism. If it doesn't know how to copy the type, it
// will give a compilation error.

#ifndef DCA_FUNCTION_COPY_FROM_HPP
#define DCA_FUNCTION_COPY_FROM_HPP

#include <cstring>
#include <complex>
#include <vector>

namespace dca {
namespace func {
// dca::func::

template <typename whatever_t>
struct copy_from {};

template <>
struct copy_from<short> {
  static void execute(int size, short* whatever_l, short* whatever_r) {
    std::memcpy(whatever_l, whatever_r, sizeof(short) * size);
  }
};

template <>
struct copy_from<int> {
  static void execute(int size, int* whatever_l, int* whatever_r) {
    std::memcpy(whatever_l, whatever_r, sizeof(int) * size);
  }
};

template <>
struct copy_from<long> {
  static void execute(int size, long* whatever_l, long* whatever_r) {
    std::memcpy(whatever_l, whatever_r, sizeof(long) * size);
  }
};

template <>
struct copy_from<std::size_t> {
  static void execute(int size, std::size_t* whatever_l, std::size_t* whatever_r) {
    std::memcpy(whatever_l, whatever_r, sizeof(std::size_t) * size);
  }
};

template <>
struct copy_from<float> {
  static void execute(int size, float* whatever_l, float* whatever_r) {
    std::memcpy(whatever_l, whatever_r, sizeof(float) * size);
  }
};

template <>
struct copy_from<double> {
  static void execute(int size, double* whatever_l, double* whatever_r) {
    std::memcpy(whatever_l, whatever_r, sizeof(double) * size);
  }
};

template <>
struct copy_from<std::complex<float>> {
  static void execute(int size, std::complex<float>* whatever_l, std::complex<float>* whatever_r) {
    std::memcpy(whatever_l, whatever_r, sizeof(std::complex<float>) * size);
  }
};

template <>
struct copy_from<std::complex<double>> {
  static void execute(int size, std::complex<double>* whatever_l, std::complex<double>* whatever_r) {
    std::memcpy(whatever_l, whatever_r, sizeof(std::complex<double>) * size);
  }
};

template <typename whatever_t>
struct copy_from<std::vector<whatever_t>> {
  static void execute(int size, std::vector<whatever_t>* whatever_l,
                      std::vector<whatever_t>* whatever_r) {
    for (int l = 0; l < size; l++) {
      whatever_l[l] = whatever_r[l];
    }
  }
};

}  // func
}  // dca

#endif  // DCA_FUNCTION_COPY_FROM_HPP
