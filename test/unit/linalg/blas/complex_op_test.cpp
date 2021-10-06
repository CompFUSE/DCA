// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//

#include <complex>
#include "gtest/gtest.h"
#include "complex_op_test_kernels.hpp"

class ComplexOpGPUTest: public::testing::Test {
public:
  double blah;
};

typedef ::testing::Types<std::complex<float>, std::complex<double>> ComplexTypes;

TYPED_TEST_CASE(ComplexOpGPUTest, ComplexTypes);

TYPED_TEST(ComplexOpGPUTest, opmult) {
  using Scalar = TypeParam;
  Scalar a{1,1};
  Scalar b{1, -2};
  gpu_oprator_opmult(a, b);
}
