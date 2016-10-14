// Copyright (C)-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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

#ifdef DCA_HAVE_CUDA
// CUBLAS
#include "dca/linalg/blas/cublas1.hpp"
#include "dca/linalg/blas/cublas3.hpp"
#endif  // DCA_HAVE_CUDA

#include "dca/linalg/lapack/lapack.hpp"

// Device selector struct
#include "dca/linalg/blas/use_device.hpp"
#include "dca/linalg/lapack/use_device.hpp"
