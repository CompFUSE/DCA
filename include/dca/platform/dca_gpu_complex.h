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
 *  This file provides vender independent complex gpu headers.
 */
#ifndef DCA_GPU_COMPLEX_H
#define DCA_GPU_COMPLEX_H

#if defined(DCA_HAVE_CUDA)
#include <cuComplex.h>
#elif defined(DCA_HAVE_HIP)
#include "dca/util/cuda2hip.h"
#include <hip/hip_complex.h>
#endif

#endif
