// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides vender independent basic gpu headers.
 *
 *  This file turns our to be essential to make sure compilation units in libraries
 *  include the same haves_defines.hpp and have the expected symbols defined.
 *
 *  Since having DCA_HAVE_GPU defined means at least basic GPU types need to be defined
 *  this is more often included rather than haves_defines.hpp directly
 */
#ifndef DCA_GPU_H
#define DCA_GPU_H

#include "dca/config/haves_defines.hpp"
#if defined(DCA_HAVE_CUDA)
#include <cuda.h>
#include <cuda_runtime.h>
#include "dca/linalg/util/error_cuda.hpp"
#elif defined(DCA_HAVE_HIP)
#include <hip/hip_runtime.h>
#include "dca/util/cuda2hip.h"
#include "dca/linalg/util/error_hip.hpp"
#endif

#endif
