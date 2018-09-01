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
  EXPECT_LE(1, stream_cont.get_max_threads());

  stream_cont.resize(7);
  EXPECT_EQ(7, stream_cont.get_max_threads());

  const int max_threads = 4;
  // initializeStreamContainer resizes the static getStreamContainer() object.
  dca::linalg::util::initializeStreamContainer(max_threads);
  EXPECT_EQ(max_threads, stream_cont.get_max_threads());

  const int streams_per_thread = stream_cont.get_streams_per_thread();
  for (int t_id = 0; t_id < max_threads; ++t_id) {
    for (int s_id = 0; s_id < streams_per_thread; ++s_id) {
      EXPECT_EQ(stream_cont(t_id, s_id), dca::linalg::util::getStream(t_id, s_id));
      EXPECT_FALSE(NULL == stream_cont(t_id, s_id));
    }
  }
}

TEST(StreamHandelsTest, Handels) {
  auto& handle_cont = dca::linalg::util::getHandleContainer();
  EXPECT_EQ(1, handle_cont.size());

  const int max_threads = 7;
  dca::linalg::util::initializeHandleContainer(max_threads);

  EXPECT_EQ(max_threads, handle_cont.size());
  EXPECT_LE(max_threads, dca::linalg::util::getStreamContainer().get_max_threads());

  cudaStream_t stream;
  for (int t_id = 0; t_id < max_threads; ++t_id) {
    cublasHandle_t handle = dca::linalg::util::getHandle(t_id);
    EXPECT_EQ(handle_cont[t_id], handle);
    cublasGetStream(handle, &stream);
    EXPECT_EQ(dca::linalg::util::getStream(t_id, 0), stream);

    const int streams_per_thread = dca::linalg::util::getStreamContainer().get_streams_per_thread();

    for (int s_id = 0; s_id < streams_per_thread; ++s_id) {
      cublasHandle_t handle = dca::linalg::util::getHandle(t_id, s_id);
      EXPECT_EQ(handle_cont[t_id], handle);
      cublasGetStream(handle, &stream);
      EXPECT_EQ(dca::linalg::util::getStream(t_id, s_id), stream);
    }
  }
}
