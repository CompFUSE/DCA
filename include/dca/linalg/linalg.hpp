// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file includes all the header files in include/dca/linalg.
// TODO: This file is temporary and will be removed or updated.

#include "dca/linalg/vector.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"

// BLAS
#include "dca/linalg/blas/blas1.hpp"
#include "dca/linalg/blas/blas2.hpp"
#include "dca/linalg/blas/blas3.hpp"

#include "dca/linalg/lapack/bennet_update.hpp"
#include "dca/linalg/lapack/inverse.hpp"
#include "dca/linalg/lapack/lapack.hpp"
#include "dca/linalg/lapack/solve.hpp"

#ifdef DCA_HAVE_CUDA
// CUBLAS
#include "dca/linalg/blas/cublas1.hpp"
#include "dca/linalg/blas/cublas3.hpp"
#include "dca/linalg/blas/cublas_conversion_char_types.hpp"
#include "dca/linalg/blas/kernels_gpu.hpp"

#include "dca/linalg/lapack/laset_gpu.hpp"
#include "dca/linalg/lapack/magma.hpp"
#include "dca/linalg/lapack/multiply_diagonal_gpu.hpp"
#endif  // DCA_HAVE_CUDA

// Device selector struct
#include "dca/linalg/device_type.hpp"
#include "dca/linalg/blas/use_device.hpp"
#include "dca/linalg/lapack/use_device.hpp"

// Utils
#include "dca/linalg/util/allocators/allocators.hpp"
#include "dca/linalg/util/copy.hpp"
#include "dca/linalg/util/lapack_exception.hpp"
#include "dca/linalg/util/util_lapack.hpp"
#include "dca/linalg/util/util_matrixop.hpp"
#ifdef DCA_HAVE_CUDA
#include "dca/linalg/util/cast_cuda.hpp"
#include "dca/linalg/util/error_cublas.hpp"
#include "dca/linalg/util/error_cuda.hpp"
#include "dca/linalg/util/handle_functions.hpp"
#include "dca/linalg/util/info_cuda.hpp"
#include "dca/linalg/util/stream_container.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/linalg/util/util_cublas.hpp"
#endif  // DCA_HAVE_CUDA
