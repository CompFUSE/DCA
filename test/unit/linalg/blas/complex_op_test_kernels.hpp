#ifndef DCA_COMPLEX_OP_TEST_KERNELS_HPP
#define DCA_COMPLEX_OP_TEST_KERNELS_HPP

#include <complex>
#include "dca/config/haves_defines.hpp"
#include "dca/platform/dca_gpu.h"
#include "dca/platform/dca_gpu_complex.h"
#include "dca/linalg/util/cast_gpu.hpp"

template <typename Type>
void gpu_operator_opmult(Type* a, const Type* b);

template <typename Type>
void gpu_operator_opmult(std::complex<Type>* a, const std::complex<Type>* b)
{
  auto* cu_a = castCudaComplex(a);
  auto* cu_b = castCudaComplex(b);
  gpu_operator_optmult(a,b);
}

#endif
