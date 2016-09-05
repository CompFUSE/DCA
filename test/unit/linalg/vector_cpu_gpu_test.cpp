// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the Vector<CPU> class.

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
  int thread_id = 2;
  int stream_id = 5;

  // Tests all the constructors.
  {
    dca::linalg::Vector<float, dca::linalg::CPU> vec(name, size, capacity, thread_id, stream_id);
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
  {
    dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> vec(name, size, capacity);
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
  ASSERT_EQ(vec.get_name(), vec_copy.get_name());
  ASSERT_EQ(vec.size(), vec_copy.size());
  ASSERT_LE(vec.size(), vec_copy.capacity());
  ASSERT_TRUE(testing::isDevicePointer(vec_copy.ptr()));

  dca::linalg::Vector<float, dca::linalg::CPU> vec_copy_copy(vec_copy);
  EXPECT_EQ(vec.get_name(), vec_copy_copy.get_name());
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
    ASSERT_EQ(vec.get_name(), vec_copy.get_name());
    ASSERT_EQ(vec.size(), vec_copy.size());
    ASSERT_EQ(capacity, vec_copy.capacity());
    ASSERT_EQ(old_ptr, vec_copy.ptr());
    ASSERT_TRUE(testing::isDevicePointer(vec_copy.ptr()));

    vec_copy_copy = vec_copy;
    EXPECT_EQ(vec.get_name(), vec_copy_copy.get_name());
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
    ASSERT_EQ(vec.get_name(), vec_copy.get_name());
    ASSERT_EQ(vec.size(), vec_copy.size());
    ASSERT_LE(vec.size(), vec_copy.capacity());
    ASSERT_FALSE(testing::isHostPointer(vec_copy.ptr()));

    vec_copy_copy = vec_copy;
    EXPECT_EQ(vec.get_name(), vec_copy_copy.get_name());
    EXPECT_EQ(vec.size(), vec_copy_copy.size());
    EXPECT_LE(vec.size(), vec_copy_copy.capacity());
    EXPECT_TRUE(testing::isHostPointer(vec_copy_copy.ptr()));

    for (int i = 0; i < vec.size(); ++i) {
      EXPECT_EQ(vec[i], vec_copy_copy[i]);
      EXPECT_NE(vec.ptr(i), vec_copy_copy.ptr(i));
    }
  }
}
