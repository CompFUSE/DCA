// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the handles and streams related api.

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
  // resizeStreamContainer resizes the global getStreamContainer() object.
  dca::linalg::util::resizeStreamContainer(max_threads);
  EXPECT_EQ(max_threads, stream_cont.get_max_threads());

  const int streams_per_thread = stream_cont.get_streams_per_thread();
  for (int t_id = 0; t_id < max_threads; ++t_id) {
    for (int s_id = 0; s_id < streams_per_thread; ++s_id) {
      EXPECT_EQ(stream_cont(t_id, s_id), dca::linalg::util::getStream(t_id, s_id));
      EXPECT_FALSE(NULL == stream_cont(t_id, s_id));
    }
  }
}

TEST(StreamHandelsTest, Handles) {
  auto& handle_cont = dca::linalg::util::getHandleContainer();
  EXPECT_EQ(1, handle_cont.size());

  auto check_stream = [&handle_cont](const int thread_id) {
    cudaStream_t stream;
    cublasGetStream(handle_cont[thread_id], &stream);

    const int streams_per_thread = dca::linalg::util::getStreamContainer().get_streams_per_thread();
    for (int s_id = 0; s_id < streams_per_thread; ++s_id) {
      cublasHandle_t handle = dca::linalg::util::getHandle(thread_id, s_id);
      EXPECT_EQ(handle_cont[thread_id], handle);
      cublasGetStream(handle, &stream);
      EXPECT_EQ(dca::linalg::util::getStream(thread_id, s_id), stream);
    }
  };

  check_stream(0);

  const int max_threads = 7;
  dca::linalg::util::resizeHandleContainer(max_threads);

  EXPECT_EQ(max_threads, handle_cont.size());
  EXPECT_LE(max_threads, dca::linalg::util::getStreamContainer().get_max_threads());

  for (int t_id = 0; t_id < max_threads; ++t_id)
    check_stream(t_id);
}
