// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file pulls in all Linalg *.h files.
// It is self-contained and can be included whenever a file depends on any of the Linalg *.h files.
//
// TODO: Make all header files self-contained.

#ifndef COMP_LIBRARY_LINALG_LINALG_HPP
#define COMP_LIBRARY_LINALG_LINALG_HPP

#include <cassert>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/linalg/linalg.hpp"

#include "comp_library/blas_lapack_plans/blas_lapack_plans.hpp"

#include "comp_library/linalg/linalg_device_types.h"

// LAPACK
#include "src/linalg_operations/GESV_tem.h"
#include "src/linalg_operations/GESV_CPU.h"

// Performance inspector
#include "src/linalg_structures/performance_inspector_tem.h"
#include "src/linalg_structures/performance_inspector_CPU.h"
#ifdef DCA_HAVE_CUDA
#include "src/linalg_structures/performance_inspector_GPU.h"
#endif

#endif  // COMP_LIBRARY_LINALG_LINALG_HPP
