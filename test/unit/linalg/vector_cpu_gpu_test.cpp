// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the interaction between Vector<CPU> and Vector<GPU>.

#include "dca/linalg/vector.hpp"
#include <complex>
#include <string>
#include <utility>
#include "gtest/gtest.h"
#include "gpu_test_util.hpp"

TEST(VectorCPUTest, PointerMemoryType) {
  size_t size = 3;
  size_t capacity = 11;
  std::string name("vector name");

  // Tests all the constructors.
  {
    dca::linalg::Vector<float, dca::linalg::CPU> vec(name, size, capacity);
    ASSERT_TRUE(testing::isHostPointer(vec.ptr()));
  }
  {
    dca::linalg::Vector<int, dca::linalg::CPU> vec(size);
    EXPECT_TRUE(testing::isHostPointer(vec.ptr()));
  }
  {
    dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> vec(size, capacity);
    EXPECT_TRUE(testing::isHostPointer(vec.ptr()));
  }
  {
    dca::linalg::Vector<int, dca::linalg::CPU> vec(name, size);
    EXPECT_TRUE(testing::isHostPointer(vec.ptr()));
  }
}

TEST(VectorCPUGPUTest, Constructors) {
  size_t size = 3;

  dca::linalg::Vector<float, dca::linalg::CPU> vec("name", size);
  // Set the elements.
  for (int i = 0; i < vec.size(); ++i) {
    float el = 3 * i - 2;
    vec[i] = el;
  }

  dca::linalg::Vector<float, dca::linalg::GPU> vec_copy(vec);
  ASSERT_EQ(vec.size(), vec_copy.size());
  ASSERT_LE(vec.size(), vec_copy.capacity());
  ASSERT_TRUE(testing::isDevicePointer(vec_copy.ptr()));

  dca::linalg::Vector<float, dca::linalg::CPU> vec_copy_copy(vec_copy);
  EXPECT_EQ(vec.size(), vec_copy_copy.size());
  EXPECT_LE(vec.size(), vec_copy_copy.capacity());
  EXPECT_TRUE(testing::isHostPointer(vec_copy_copy.ptr()));

  for (int i = 0; i < vec.size(); ++i) {
    EXPECT_EQ(vec[i], vec_copy_copy[i]);
    EXPECT_NE(vec.ptr(i), vec_copy_copy.ptr(i));
  }
}

TEST(VectorCPUGPUTest, Assignement) {
  {
    // Assign a vector that fits into the capacity.
    size_t size = 3;

    dca::linalg::Vector<float, dca::linalg::GPU> vec_copy(10);
    auto old_ptr = vec_copy.ptr();
    auto capacity = vec_copy.capacity();
    dca::linalg::Vector<float, dca::linalg::CPU> vec_copy_copy(6);
    auto old_ptr_2 = vec_copy_copy.ptr();
    auto capacity_2 = vec_copy_copy.capacity();

    dca::linalg::Vector<float, dca::linalg::CPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy = vec;
    ASSERT_EQ(vec.size(), vec_copy.size());
    ASSERT_EQ(capacity, vec_copy.capacity());
    ASSERT_EQ(old_ptr, vec_copy.ptr());
    ASSERT_TRUE(testing::isDevicePointer(vec_copy.ptr()));

    vec_copy_copy = vec_copy;
    EXPECT_EQ(vec.size(), vec_copy_copy.size());
    EXPECT_EQ(capacity_2, vec_copy_copy.capacity());
    EXPECT_EQ(old_ptr_2, vec_copy_copy.ptr());
    EXPECT_TRUE(testing::isHostPointer(vec_copy_copy.ptr()));

    for (int i = 0; i < vec.size(); ++i) {
      EXPECT_EQ(vec[i], vec_copy_copy[i]);
      EXPECT_NE(vec.ptr(i), vec_copy_copy.ptr(i));
    }
  }
  {
    // Assign a vector that doesn't fit into the capacity.
    dca::linalg::Vector<float, dca::linalg::GPU> vec_copy(10);
    dca::linalg::Vector<float, dca::linalg::CPU> vec_copy_copy(6);
    size_t size = std::max(vec_copy.capacity(), vec_copy_copy.capacity()) + 1;

    dca::linalg::Vector<float, dca::linalg::CPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy = vec;
    ASSERT_EQ(vec.size(), vec_copy.size());
    ASSERT_LE(vec.size(), vec_copy.capacity());
    ASSERT_FALSE(testing::isHostPointer(vec_copy.ptr()));

    vec_copy_copy = vec_copy;
    EXPECT_EQ(vec.size(), vec_copy_copy.size());
    EXPECT_LE(vec.size(), vec_copy_copy.capacity());
    EXPECT_TRUE(testing::isHostPointer(vec_copy_copy.ptr()));

    for (int i = 0; i < vec.size(); ++i) {
      EXPECT_EQ(vec[i], vec_copy_copy[i]);
      EXPECT_NE(vec.ptr(i), vec_copy_copy.ptr(i));
    }
  }
}

TEST(VectorCPUGPUTest, Set) {
  {
    // Assign a vector that fits into the capacity.
    size_t size = 3;

    dca::linalg::Vector<float, dca::linalg::GPU> vec_copy(10);
    auto old_ptr = vec_copy.ptr();
    auto capacity = vec_copy.capacity();
    dca::linalg::Vector<float, dca::linalg::CPU> vec_copy_copy(6);
    auto old_ptr_2 = vec_copy_copy.ptr();
    auto capacity_2 = vec_copy_copy.capacity();

    dca::linalg::Vector<float, dca::linalg::CPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy.set(vec, 0, 1);
    ASSERT_EQ(vec.size(), vec_copy.size());
    ASSERT_EQ(capacity, vec_copy.capacity());
    ASSERT_EQ(old_ptr, vec_copy.ptr());
    ASSERT_TRUE(testing::isDevicePointer(vec_copy.ptr()));

    vec_copy_copy.set(vec_copy, 0, 1);
    EXPECT_EQ(vec.size(), vec_copy_copy.size());
    EXPECT_EQ(capacity_2, vec_copy_copy.capacity());
    EXPECT_EQ(old_ptr_2, vec_copy_copy.ptr());
    EXPECT_TRUE(testing::isHostPointer(vec_copy_copy.ptr()));

    for (int i = 0; i < vec.size(); ++i) {
      EXPECT_EQ(vec[i], vec_copy_copy[i]);
      EXPECT_NE(vec.ptr(i), vec_copy_copy.ptr(i));
    }
  }
  {
    // Assign a vector that doesn't fit into the capacity.
    dca::linalg::Vector<float, dca::linalg::GPU> vec_copy(10);
    dca::linalg::Vector<float, dca::linalg::CPU> vec_copy_copy(6);
    size_t size = std::max(vec_copy.capacity(), vec_copy_copy.capacity()) + 1;

    dca::linalg::Vector<float, dca::linalg::CPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy.set(vec, 0, 1);
    ASSERT_EQ(vec.size(), vec_copy.size());
    ASSERT_LE(vec.size(), vec_copy.capacity());
    ASSERT_FALSE(testing::isHostPointer(vec_copy.ptr()));

    vec_copy_copy.set(vec_copy, 0, 1);
    EXPECT_EQ(vec.size(), vec_copy_copy.size());
    EXPECT_LE(vec.size(), vec_copy_copy.capacity());
    EXPECT_TRUE(testing::isHostPointer(vec_copy_copy.ptr()));

    for (int i = 0; i < vec.size(); ++i) {
      EXPECT_EQ(vec[i], vec_copy_copy[i]);
      EXPECT_NE(vec.ptr(i), vec_copy_copy.ptr(i));
    }
  }
}

TEST(VectorCPUTest, setAsync) {
  std::vector<int> vec(4, 1);

  dca::linalg::Vector<int, dca::linalg::GPU> vec_copy;
  dca::linalg::Vector<int, dca::linalg::CPU> vec_copy_copy;

  cudaStream_t stream;
  cudaStreamCreate(&stream);

  vec_copy.setAsync(vec, stream);
  vec_copy_copy.setAsync(vec_copy, stream);
  cudaStreamSynchronize(stream);

  EXPECT_EQ(vec.size(), vec_copy_copy.size());
  for (int i = 0; i < vec.size(); ++i)
    EXPECT_EQ(vec[i], vec_copy_copy[i]);

  cudaStreamDestroy(stream);
}
