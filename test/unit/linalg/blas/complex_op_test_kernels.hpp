// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file tests whether the platform has working assignment operators *= etc.

#ifndef DCA_COMPLEX_OP_TEST_KERNELS_HPP
#define DCA_COMPLEX_OP_TEST_KERNELS_HPP

#include <complex>
#include "dca/config/haves_defines.hpp"
#include "dca/platform/dca_gpu.h"
#include "dca/platform/dca_gpu_complex.h"
#include "dca/linalg/util/cast_gpu.hpp"
#include "dca/linalg/util/stream_functions.hpp"
template <typename Type>
void gpu_operator_opmult(Type* a, const Type* b, int thread_id, int stream_id);

extern template void gpu_operator_opmult(magmaFloatComplex* a, const magmaFloatComplex* b, int thread_id , int stream_id);
extern template void gpu_operator_opmult(magmaDoubleComplex* a, const magmaDoubleComplex* b, int thread_id , int stream_id);

template <typename Type>
void gpu_operator_opmult(std::complex<Type>* a, const std::complex<Type>* b, int thread_id, int stream_id)
{
  using namespace dca::linalg::util;
  auto* cu_a = castCudaComplex(a);
  auto* cu_b = castCudaComplex(b);
  gpu_operator_opmult(cu_a,cu_b, thread_id, stream_id);
}


#endif
