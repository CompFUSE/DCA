// Copyright (C)-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file pulls in all BLAS and LAPACK plans *.h files.
// It is self-contained and can be included whenever a file depends on any of the BLAS and LAPACK
// plans *.h files.
//
// TODO: Make all header files self-contained.

// BLAS
#include "dca/linalg/blas/blas1.hpp"
#include "dca/linalg/blas/blas2.hpp"
#include "dca/linalg/blas/blas3.hpp"

#ifdef DCA_HAVE_CUDA
// CUBLAS
#include "dca/linalg/blas/cublas1.hpp"
#include "dca/linalg/blas/cublas3.hpp"
#endif  // DCA_HAVE_CUDA

// Device selector struct
#include "dca/linalg/blas/usedevice.hpp"
