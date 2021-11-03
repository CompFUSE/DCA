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
#include "dca/linalg/vector.hpp"
#include "complex_op_test_kernels.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/platform/dca_gpu.h"
#include "dca/linalg/util/gpu_stream.hpp"

template<typename Type>
class ComplexOpGPUTest: public::testing::Test {
public:
  double blah;
};

typedef ::testing::Types<std::complex<float>, std::complex<double>> ComplexTypes;

TYPED_TEST_CASE(ComplexOpGPUTest, ComplexTypes);

TYPED_TEST(ComplexOpGPUTest, opmult) {
  using dca::linalg::DeviceType;
  using Scalar = TypeParam;
  dca::linalg::Vector<Scalar, DeviceType::GPU> gpu_vec(2);
  dca::linalg::Vector<Scalar, DeviceType::CPU> cpu_vec(2);
  Scalar a{1,1};
  Scalar b{1, -2};
  cpu_vec[0] = {1,1};
  cpu_vec[1] = {1, -2};
  gpu_vec.setAsync(cpu_vec, 0);
  auto dev_ptr = gpu_vec.ptr();
  
  gpu_operator_opmult_cast(dev_ptr, dev_ptr+1, 0, 0);
}
