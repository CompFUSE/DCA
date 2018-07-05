// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Redefine keyword like __device__ and __host__
// so that code can be compiled either by nvcc or the host compiler.

#ifndef DCA_UTIL_CUDA_DEFINITIONS_HPP
#define DCA_UTIL_CUDA_DEFINITIONS_HPP

#ifdef __CUDACC__
#define __HOST__ __host__
#define __DEVICE__ __device__
#define IS_CUDA_COMPILER

#else // the compiler is not nvcc
#define __HOST__
#define __DEVICE__
#endif


#endif  // DCA_UTIL_CUDA_DEFINITIONS_HPP
