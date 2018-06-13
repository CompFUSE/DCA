// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the return code checking of Cuda and Cublas.

#include "dca/linalg/util/error_cublas.hpp"
#include "dca/linalg/util/error_cuda.hpp"
#include <stdexcept>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "gtest/gtest.h"

TEST(CheckRCTest, Cuda) {
  cudaError_t codes[] = {cudaErrorDeviceAlreadyInUse, cudaErrorIllegalAddress,
                         cudaErrorIllegalInstruction, cudaErrorInvalidDevice,
                         cudaErrorInvalidPitchValue,  cudaErrorInvalidValue,
                         cudaErrorLaunchFailure,      cudaErrorMemoryAllocation,
                         cudaErrorNoDevice,           cudaErrorUnknown};

  EXPECT_NO_THROW(checkRC(cudaSuccess));

  for (cudaError_t err : codes) {
    EXPECT_THROW(checkRC(err), std::logic_error);
  }
}

TEST(CheckRCTest, Cublas) {
  cublasStatus_t codes[] = {CUBLAS_STATUS_NOT_INITIALIZED, CUBLAS_STATUS_INVALID_VALUE,
                            CUBLAS_STATUS_EXECUTION_FAILED, CUBLAS_STATUS_NOT_SUPPORTED};

  EXPECT_NO_THROW(checkRC(CUBLAS_STATUS_SUCCESS));

  for (cublasStatus_t err : codes) {
    EXPECT_THROW(checkRC(err), std::logic_error);
  }
}
