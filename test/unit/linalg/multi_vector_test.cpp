// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (galduzz@itp.phys.ethz.ch)
//
// This file tests the MultiVector.

#include "dca/linalg/multi_vector.hpp"

#include <vector>

#include "gtest/gtest.h"

TEST(MultiVector, ResizeAndGet) {
  dca::linalg::MultiVector<dca::linalg::CPU, int, double, char> v;
  v.resizeNoCopy(2);

  const std::vector<int> v_int{0, 1};
  const std::vector<double> v_double{0.5, 1.5};
  const std::vector<char> v_char{'a', 'b'};

  for (int i = 0; i < v.size(); ++i) {
    v.get<0>()[i] = v_int[i];
    v.get<1>()[i] = v_double[i];
    v.get<2>()[i] = v_char[i];
  }

  for (int i = 0; i < v.size(); ++i) {
    EXPECT_EQ(v.get<0>()[i], v_int[i]);
    EXPECT_EQ(v.get<1>()[i], v_double[i]);
    EXPECT_EQ(v.get<2>()[i], v_char[i]);
  }
}
