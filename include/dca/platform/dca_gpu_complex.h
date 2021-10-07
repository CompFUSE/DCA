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
 *  This file provides working vender or magma complex gpu headers.
 */
#ifndef DCA_GPU_COMPLEX_H
#define DCA_GPU_COMPLEX_H

#if defined(DCA_HAVE_CUDA)
#include <cuComplex.h>
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#elif defined(DCA_HAVE_HIP)
// hipComplex types are faulty so we use the magma complex types and operators
#include <magma_operators.h>
#include "dca/util/cuda2hip.h"
#endif

#endif
