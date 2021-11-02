// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//

/** \file
 *  This file tests the return code checking of Hip or Cuda and HipBLAS or CuBLAS.
 */
#include "dca/platform/dca_gpu.h"
#include "dca/platform/dca_gpu_blas.h"
#include <stdexcept>
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
