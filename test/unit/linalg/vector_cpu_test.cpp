// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the Vector<CPU> class.

#include "dca/linalg/vector.hpp"
#include <complex>
#include <string>
#include <utility>
#include "gtest/gtest.h"

TEST(VectorCPUTest, Constructors) {
  size_t size = 3;
  size_t capacity = 11;
  std::string name("vector name");

  // Tests all the constructors.
  {
    dca::linalg::Vector<float, dca::linalg::CPU> vec(name, size, capacity);
    ASSERT_EQ(name, vec.get_name());
    ASSERT_EQ(size, vec.size());
    ASSERT_LE(capacity, vec.capacity());
    ASSERT_NE(nullptr, vec.ptr());
  }
  {
    dca::linalg::Vector<double, dca::linalg::CPU> vec;
    EXPECT_EQ(0, vec.size());
    EXPECT_LE(0, vec.capacity());
  }
  {
    dca::linalg::Vector<int, dca::linalg::CPU> vec(size);
    EXPECT_EQ(size, vec.size());
    EXPECT_LE(size, vec.capacity());
    EXPECT_NE(nullptr, vec.ptr());
  }
  {
    dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> vec(size, capacity);
    EXPECT_EQ(size, vec.size());
    EXPECT_LE(capacity, vec.capacity());
    EXPECT_NE(nullptr, vec.ptr());
  }
  {
    dca::linalg::Vector<double, dca::linalg::CPU> vec(name);
    EXPECT_EQ(name, vec.get_name());
    EXPECT_EQ(0, vec.size());
    EXPECT_LE(0, vec.capacity());
  }
  {
    dca::linalg::Vector<int, dca::linalg::CPU> vec(name, size);
    EXPECT_EQ(name, vec.get_name());
    EXPECT_EQ(size, vec.size());
    EXPECT_LE(size, vec.capacity());
    EXPECT_NE(nullptr, vec.ptr());
  }
}

TEST(VectorCPUTest, ElementPointers) {
  // Check if the pointers are computed correctly.
  size_t size = 5;

  dca::linalg::Vector<int, dca::linalg::CPU> vec(size);
  const dca::linalg::Vector<int, dca::linalg::CPU>& vec_const_ref(vec);
  for (int i = 0; i < vec.size(); ++i) {
    int* ptr = vec.ptr();
    EXPECT_EQ(i, vec.ptr(i) - ptr);
    EXPECT_EQ(vec.ptr(i), vec_const_ref.ptr(i));
    EXPECT_EQ(vec.ptr(i), &vec[i]);
    EXPECT_EQ(vec.ptr(i), &vec_const_ref[i]);
  }
}

TEST(VectorCPUTest, ElementAccess) {
  // Check if the different element accesses return the same value.
  size_t size = 4;

  dca::linalg::Vector<int, dca::linalg::CPU> vec(size);
  const dca::linalg::Vector<int, dca::linalg::CPU>& vec_const_ref(vec);
  for (int i = 0; i < vec.size(); ++i) {
    int el = 3 * i - 2;
    vec[i] = el;
    EXPECT_EQ(el, vec[i]);
    EXPECT_EQ(el, vec_const_ref[i]);
    EXPECT_EQ(el, *(vec.ptr(i)));
    EXPECT_EQ(el, *(vec_const_ref.ptr(i)));
  }
}

TEST(VectorCPUTest, CopyConstructor) {
  size_t size = 4;

  dca::linalg::Vector<float, dca::linalg::CPU> vec("name", size);
  // Set the elements.
  for (int i = 0; i < vec.size(); ++i) {
    float el = 3 * i - 2;
    vec[i] = el;
  }

  dca::linalg::Vector<float, dca::linalg::CPU> vec_copy(vec, "another name");
  EXPECT_EQ("another name", vec_copy.get_name());
  EXPECT_EQ(vec.size(), vec_copy.size());
  EXPECT_LE(vec.size(), vec_copy.capacity());

  for (int i = 0; i < vec.size(); ++i) {
    EXPECT_EQ(vec[i], vec_copy[i]);
    EXPECT_NE(vec.ptr(i), vec_copy.ptr(i));
  }
}

TEST(VectorCPUTest, Assignement) {
  {
    // Assign a vector that fits into the capacity.
    size_t size = 4;

    dca::linalg::Vector<float, dca::linalg::CPU> vec_copy(10);
    auto old_ptr = vec_copy.ptr();
    auto capacity = vec_copy.capacity();

    dca::linalg::Vector<float, dca::linalg::CPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy = vec;
    EXPECT_EQ(vec.size(), vec_copy.size());
    EXPECT_EQ(capacity, vec_copy.capacity());
    EXPECT_EQ(old_ptr, vec_copy.ptr());

    for (int i = 0; i < vec.size(); ++i) {
      EXPECT_EQ(vec[i], vec_copy[i]);
      EXPECT_NE(vec.ptr(i), vec_copy.ptr(i));
    }
  }
  {
    // Assign a vector that does not fit into the capacity.
    dca::linalg::Vector<float, dca::linalg::CPU> vec_copy(10);
    auto size = vec_copy.capacity();
    ++size;

    dca::linalg::Vector<float, dca::linalg::CPU> vec("name", size);

    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy = vec;
    EXPECT_EQ(vec.size(), vec_copy.size());
    EXPECT_LE(vec.size(), vec_copy.capacity());

    for (int i = 0; i < vec.size(); ++i) {
      EXPECT_EQ(vec[i], vec_copy[i]);
      EXPECT_NE(vec.ptr(i), vec_copy.ptr(i));
    }
  }
}

TEST(VectorCPUTest, Set) {
  {
    // Assign a vector that fits into the capacity.
    size_t size = 4;

    dca::linalg::Vector<float, dca::linalg::CPU> vec_copy(10);
    auto old_ptr = vec_copy.ptr();
    auto capacity = vec_copy.capacity();

    dca::linalg::Vector<float, dca::linalg::CPU> vec("name", size);
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy.set(vec, 0, 1);
    EXPECT_EQ(vec.size(), vec_copy.size());
    EXPECT_EQ(capacity, vec_copy.capacity());
    EXPECT_EQ(old_ptr, vec_copy.ptr());

    for (int i = 0; i < vec.size(); ++i) {
      EXPECT_EQ(vec[i], vec_copy[i]);
      EXPECT_NE(vec.ptr(i), vec_copy.ptr(i));
    }
  }
  {
    // Assign a vector that does not fit into the capacity.
    dca::linalg::Vector<float, dca::linalg::CPU> vec_copy(10);
    auto size = vec_copy.capacity();
    ++size;

    dca::linalg::Vector<float, dca::linalg::CPU> vec("name", size);

    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      float el = 3 * i - 2;
      vec[i] = el;
    }

    vec_copy.set(vec, 0, 1);
    EXPECT_EQ(vec.size(), vec_copy.size());
    EXPECT_LE(vec.size(), vec_copy.capacity());

    for (int i = 0; i < vec.size(); ++i) {
      EXPECT_EQ(vec[i], vec_copy[i]);
      EXPECT_NE(vec.ptr(i), vec_copy.ptr(i));
    }
  }
}

TEST(VectorCPUTest, Resize) {
  {
    size_t size = 4;

    dca::linalg::Vector<long, dca::linalg::CPU> vec(size);

    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      long el = 1 + 3 * i;
      vec[i] = el;
    }

    // Resize to capacity. No reallocation has to take place.
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();
    int new_size = capacity;
    vec.resize(new_size);
    EXPECT_EQ(new_size, vec.size());
    EXPECT_EQ(capacity, vec.capacity());
    EXPECT_EQ(old_ptr, vec.ptr());

    // Check the value of the elements.
    for (int i = 0; i < size; ++i) {
      long el = 1 + 3 * i;
      EXPECT_EQ(el, vec[i]);
    }
  }
  {
    size_t size = 5;

    dca::linalg::Vector<long, dca::linalg::CPU> vec(size);
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      long el = 1 + 3 * i;
      vec[i] = el;
    }

    // Shrink the vector. No reallocation has to take place.
    int new_size = 2;
    vec.resize(new_size);
    EXPECT_EQ(new_size, vec.size());
    EXPECT_EQ(capacity, vec.capacity());
    EXPECT_EQ(old_ptr, vec.ptr());

    // Check the value of the elements.
    for (int i = 0; i < vec.size(); ++i) {
      long el = 1 + 3 * i;
      EXPECT_EQ(el, vec[i]);
    }
  }
  {
    size_t size = 3;

    dca::linalg::Vector<long, dca::linalg::CPU> vec(size);
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();
    // Set the elements.
    for (int i = 0; i < vec.size(); ++i) {
      long el = 1 + 3 * i;
      vec[i] = el;
    }

    // New size is larger than capacity.
    // Reallocation has to take place.
    int new_size = capacity + 1;
    vec.resize(new_size);
    EXPECT_EQ(new_size, vec.size());
    EXPECT_LE(new_size, vec.capacity());
    EXPECT_NE(old_ptr, vec.ptr());

    // Check the value of the elements.
    for (int i = 0; i < size; ++i) {
      long el = 1 + 3 * i;
      EXPECT_EQ(el, vec[i]);
    }
  }
}

TEST(VectorCPUTest, ResizeNoCopy) {
  {
    size_t size = 4;

    dca::linalg::Vector<long, dca::linalg::CPU> vec(size);

    // Resize to capacity. No reallocation has to take place.
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();
    size_t new_size = capacity;
    vec.resizeNoCopy(new_size);
    EXPECT_EQ(new_size, vec.size());
    EXPECT_EQ(capacity, vec.capacity());
    EXPECT_EQ(old_ptr, vec.ptr());
  }
  {
    size_t size = 5;

    dca::linalg::Vector<long, dca::linalg::CPU> vec(size);
    auto old_ptr = vec.ptr();
    auto capacity = vec.capacity();

    // Shrink the vector. No reallocation has to take place.
    size_t new_size = 2;
    vec.resizeNoCopy(new_size);
    EXPECT_EQ(new_size, vec.size());
    EXPECT_EQ(capacity, vec.capacity());
    EXPECT_EQ(old_ptr, vec.ptr());
  }
  {
    size_t size = 3;

    dca::linalg::Vector<long, dca::linalg::CPU> vec(size);
    auto capacity = vec.capacity();

    // New size is larger than capacity.
    // Reallocation has to take place.
    size_t new_size = capacity + 1;
    vec.resizeNoCopy(new_size);
    EXPECT_EQ(new_size, vec.size());
    EXPECT_LE(new_size, vec.capacity());
  }
}

TEST(VectorCPUTest, Clear) {
  dca::linalg::Vector<double, dca::linalg::CPU> vec(42);

  EXPECT_EQ(42, vec.size());
  vec.clear();
  EXPECT_EQ(0, vec.size());
  EXPECT_EQ(0, vec.capacity());
}
