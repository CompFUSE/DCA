// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the interaction between Matrix<CPU> and Matrix<GPU>.

#include "dca/linalg/matrix.hpp"
#include <complex>
#include <string>
#include <utility>
#include "gtest/gtest.h"
#include "cpu_test_util.hpp"
#include "gpu_test_util.hpp"

TEST(MatrixCPUTest, PointerMemoryType) {
  int size = 3;
  int capacity = 11;
  std::pair<int, int> size2(4, 5);
  std::pair<int, int> capacity2(13, 17);
  std::string name("matrix name");

  // Test the pointers of the constructors.
  {
    dca::linalg::Matrix<float, dca::linalg::CPU> mat(name, size2, capacity2);
    ASSERT_TRUE(testing::isHostPointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<int, dca::linalg::CPU> mat(size);
    EXPECT_TRUE(testing::isHostPointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> mat(size, capacity);
    EXPECT_TRUE(testing::isHostPointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<std::complex<float>, dca::linalg::CPU> mat(size2);
    EXPECT_TRUE(testing::isHostPointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<float, dca::linalg::CPU> mat(size2, capacity2);
    EXPECT_TRUE(testing::isHostPointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<int, dca::linalg::CPU> mat(name, size);
    EXPECT_TRUE(testing::isHostPointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> mat(name, size, capacity);
    EXPECT_TRUE(testing::isHostPointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<std::complex<float>, dca::linalg::CPU> mat(name, size2);
    EXPECT_TRUE(testing::isHostPointer(mat.ptr()));
  }
}

TEST(MatrixCPUGPUTest, Constructors) {
  std::pair<int, int> size2(2, 3);

  dca::linalg::Matrix<float, dca::linalg::CPU> mat("name", size2);
  auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
  testing::setMatrixElements(mat, el_value);

  dca::linalg::Matrix<float, dca::linalg::GPU> mat_copy(mat);

  ASSERT_EQ(mat.size(), mat_copy.size());
  ASSERT_LE(mat.size().first, mat_copy.capacity().first);
  ASSERT_LE(mat.size().second, mat_copy.capacity().second);
  ASSERT_TRUE(testing::isDevicePointer(mat_copy.ptr()));

  dca::linalg::Matrix<float, dca::linalg::CPU> mat_copy_copy(mat_copy, "another name");
  EXPECT_EQ("another name", mat_copy_copy.get_name());
  EXPECT_EQ(mat.size(), mat_copy_copy.size());
  EXPECT_LE(mat.size().first, mat_copy_copy.capacity().first);
  EXPECT_LE(mat.size().second, mat_copy_copy.capacity().second);
  EXPECT_TRUE(testing::isHostPointer(mat_copy_copy.ptr()));

  EXPECT_EQ(mat, mat_copy_copy);
}

TEST(MatrixCPUGPUTest, Assignement) {
  {
    // Assign a matrix that fits into the capacity.
    std::pair<int, int> size2(2, 3);

    dca::linalg::Matrix<float, dca::linalg::GPU> mat_copy(10);
    auto old_ptr = mat_copy.ptr();
    auto capacity = mat_copy.capacity();
    dca::linalg::Matrix<float, dca::linalg::CPU> mat_copy_copy(6);
    auto old_ptr_2 = mat_copy_copy.ptr();
    auto capacity_2 = mat_copy_copy.capacity();

    dca::linalg::Matrix<float, dca::linalg::CPU> mat("name", size2);
    auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    mat_copy = mat;
    ASSERT_EQ(mat.size(), mat_copy.size());
    ASSERT_EQ(capacity, mat_copy.capacity());
    ASSERT_EQ(old_ptr, mat_copy.ptr());
    ASSERT_TRUE(testing::isDevicePointer(mat_copy.ptr()));

    mat_copy_copy = mat_copy;
    EXPECT_EQ(mat.size(), mat_copy_copy.size());
    EXPECT_EQ(capacity_2, mat_copy_copy.capacity());
    EXPECT_EQ(old_ptr_2, mat_copy_copy.ptr());
    EXPECT_TRUE(testing::isHostPointer(mat_copy_copy.ptr()));

    for (int j = 0; j < mat.nrCols(); ++j)
      for (int i = 0; i < mat.nrRows(); ++i) {
        EXPECT_EQ(mat(i, j), mat_copy_copy(i, j));
        EXPECT_NE(mat.ptr(i, j), mat_copy_copy.ptr(i, j));
      }
  }
  {
    // Assign a matrix that doesn't into the capacity.
    dca::linalg::Matrix<float, dca::linalg::GPU> mat_copy(10);
    dca::linalg::Matrix<float, dca::linalg::CPU> mat_copy_copy(6);
    auto size2 =
        std::make_pair(std::max(mat_copy.capacity().first, mat_copy_copy.capacity().first) + 1, 3);

    dca::linalg::Matrix<float, dca::linalg::CPU> mat("name", size2);
    auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    mat_copy = mat;
    ASSERT_EQ(mat.size(), mat_copy.size());
    ASSERT_LE(mat.size().first, mat_copy.capacity().first);
    ASSERT_LE(mat.size().second, mat_copy.capacity().second);
    ASSERT_TRUE(testing::isDevicePointer(mat_copy.ptr()));

    mat_copy_copy = mat_copy;
    EXPECT_EQ(mat.size(), mat_copy_copy.size());
    EXPECT_LE(mat.size().first, mat_copy_copy.capacity().first);
    EXPECT_LE(mat.size().second, mat_copy_copy.capacity().second);
    EXPECT_TRUE(testing::isHostPointer(mat_copy_copy.ptr()));

    for (int j = 0; j < mat.nrCols(); ++j)
      for (int i = 0; i < mat.nrRows(); ++i) {
        EXPECT_EQ(mat(i, j), mat_copy_copy(i, j));
        EXPECT_NE(mat.ptr(i, j), mat_copy_copy.ptr(i, j));
      }
  }
}

TEST(MatrixCPUGPUTest, SetAsync) {
  dca::linalg::Matrix<int, dca::linalg::CPU> mat(std::make_pair(32, 30));
  dca::linalg::Matrix<int, dca::linalg::GPU> mat_copy;
  dca::linalg::Matrix<int, dca::linalg::CPU> mat_copy_copy;

  auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
  testing::setMatrixElements(mat, el_value);

  dca::linalg::util::CudaStream stream;

  mat_copy.setAsync(mat, stream);
  mat_copy_copy.setAsync(mat_copy, stream);
  stream.sync();

  EXPECT_EQ(mat, mat_copy_copy);
}
