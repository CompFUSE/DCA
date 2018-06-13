// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the handels and streams related api.

#include "dca/linalg/util/handle_functions.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include <stdexcept>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "gtest/gtest.h"

TEST(StreamHandelsTest, Streams) {
  auto& stream_cont = dca::linalg::util::getStreamContainer();

  for (int t_id = 0; t_id < DCA_MAX_THREADS; ++t_id) {
    for (int s_id = 0; s_id < DCA_STREAMS_PER_THREAD; ++s_id) {
      EXPECT_EQ(stream_cont(t_id, s_id), dca::linalg::util::getStream(t_id, s_id));
      EXPECT_FALSE(NULL == stream_cont(t_id, s_id));
    }
  }
}

TEST(StreamHandelsTest, Handels) {
  auto& handle_cont = dca::linalg::util::getHandleContainer();
  cudaStream_t stream;

  for (int t_id = 0; t_id < DCA_MAX_THREADS; ++t_id) {
    cublasHandle_t handle = dca::linalg::util::getHandle(t_id);
    EXPECT_EQ(handle_cont(t_id), handle);
    cublasGetStream(handle, &stream);
    EXPECT_EQ(dca::linalg::util::getStream(t_id, 0), stream);

    for (int s_id = 0; s_id < DCA_STREAMS_PER_THREAD; ++s_id) {
      cublasHandle_t handle = dca::linalg::util::getHandle(t_id, s_id);
      EXPECT_EQ(handle_cont(t_id), handle);
      cublasGetStream(handle, &stream);
      EXPECT_EQ(dca::linalg::util::getStream(t_id, s_id), stream);
    }
  }
}
