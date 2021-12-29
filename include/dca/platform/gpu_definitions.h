// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//          Peter Doak (doakpw@ornl.gov)
//
// Redefine keyword like __device__ and __host__
// so that code can be compiled by hipcc, nvcc or the host compiler.

#ifndef DCA_UTIL_GPU_DEFINITIONS_HPP
#define DCA_UTIL_GPU_DEFINITIONS_HPP

#include "dca/config/haves_defines.hpp"

#ifdef DCA_HAVE_GPU
#include "dca/platform/dca_gpu.h"
#endif

// clang-format off
#if defined(DCA_HAVE_CUDA)
  #ifdef __CUDACC__
    #define __HOST__ __host__
    #define __DEVICE__ __device__
    #define __CONSTANT__ __constant__
    #define IS_GPU_COMPILER
  #else // the compiler is an cpu compiler
    #define __HOST__
    #define __DEVICE__
    #define __CONSTANT__
  #endif
#elif defined(DCA_HAVE_HIP)
  #ifdef __HIP_PLATFORM_AMD__
    #define __HOST__ __host__
    #define __DEVICE__ __device__
    #define __CONSTANT__ __constant__
  #else // the compiler is an cpu compiler
    #define __HOST__
    #define __DEVICE__
    #define __CONSTANT__
  #endif
#else // no GPU type defined
    #define __HOST__
    #define __DEVICE__
    #define __CONSTANT__
#endif
// clang-format on


#endif  // DCA_UTIL_CUDA_DEFINITIONS_HPP
