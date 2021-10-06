// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file provides the castCudaComplex utility and the CudaComplex typedef.


#include <cassert>
#include "complex_op_test_kernels.hpp"

template<typename Type>
__global__ void gpu_operator_opmult_kernel(Type* a, const Type* b) {
  a[0] *= b[0];
}

template <typename Type>
void gpu_operator_opmult(Type* a, const Type* b) {
gpu_operator_opmult_kernel<<1, 1, 0, stream>>(a, b);
}
