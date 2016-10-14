//-*-C++-*-

#include <complex>
#include <iostream>
#include <stdexcept>
#include <cstdio>

#include "linalg_device_types.h"

#include "cuda_runtime.h"
#include "cublas_v2.h"

#define HAVE_CUBLAS 1
#include "magma.h"

#include "dca/util/integer_division.hpp"
#include "dca/linalg/util/error_cuda.hpp"
#include "dca/linalg/util/error_cublas.hpp"
#include "dca/linalg/util/handle_functions.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/linalg/util/util_cublas.hpp"

#include "basic_cuda_functions.cu.h"

#include <utility>

// own kernels
#include "src/linalg_operations/GEMD_GPU.cu.h"

// magma
#include "src/linalg_operations/GETRF_GPU.cu.h"
#include "src/linalg_operations/GETRI_GPU.cu.h"