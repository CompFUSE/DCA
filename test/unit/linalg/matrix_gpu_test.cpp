// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the Matrix<GPU> class.

#include "dca/linalg/matrix.hpp"
#include <complex>
#include <string>
#include <utility>
#include "gtest/gtest.h"
#include "gpu_test_util.hpp"

TEST(MatrixGPUTest, Constructors) {
  int size = 3;
  int capacity = 11;
  std::pair<int, int> size2(4, 5);
  std::pair<int, int> capacity2(13, 17);
  std::string name("matrix name");

  // Tests all the constructors.
  {
    dca::linalg::Matrix<float, dca::linalg::GPU> mat(name, size2, capacity2);
    ASSERT_EQ(name, mat.get_name());
    ASSERT_EQ(size2, mat.size());
    ASSERT_LE(capacity2.first, mat.capacity().first);
    ASSERT_LE(capacity2.second, mat.capacity().second);
    ASSERT_NE(nullptr, mat.ptr());
    ASSERT_TRUE(testing::isDevicePointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<double, dca::linalg::GPU> mat;
    EXPECT_EQ(std::make_pair(0, 0), mat.size());
    EXPECT_LE(0, mat.capacity().first);
    EXPECT_LE(0, mat.capacity().second);
  }
  {
    dca::linalg::Matrix<int, dca::linalg::GPU> mat(size);
    EXPECT_EQ(std::make_pair(size, size), mat.size());
    EXPECT_LE(size, mat.capacity().first);
    EXPECT_LE(size, mat.capacity().second);
    EXPECT_NE(nullptr, mat.ptr());
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<std::complex<double>, dca::linalg::GPU> mat(size, capacity);
    EXPECT_EQ(std::make_pair(size, size), mat.size());
    EXPECT_LE(capacity, mat.capacity().first);
    EXPECT_LE(capacity, mat.capacity().second);
    EXPECT_NE(nullptr, mat.ptr());
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<std::complex<float>, dca::linalg::GPU> mat(size2);
    EXPECT_EQ(size2, mat.size());
    EXPECT_LE(size2.first, mat.capacity().first);
    EXPECT_LE(size2.second, mat.capacity().second);
    EXPECT_NE(nullptr, mat.ptr());
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<float, dca::linalg::GPU> mat(size2, capacity2);
    EXPECT_EQ(size2, mat.size());
    EXPECT_LE(capacity2.first, mat.capacity().first);
    EXPECT_LE(capacity2.second, mat.capacity().second);
    EXPECT_NE(nullptr, mat.ptr());
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<double, dca::linalg::GPU> mat(name);
    EXPECT_EQ(name, mat.get_name());
    EXPECT_EQ(std::make_pair(0, 0), mat.size());
    EXPECT_LE(0, mat.capacity().first);
    EXPECT_LE(0, mat.capacity().second);
  }
  {
    dca::linalg::Matrix<int, dca::linalg::GPU> mat(name, size);
    EXPECT_EQ(name, mat.get_name());
    EXPECT_EQ(std::make_pair(size, size), mat.size());
    EXPECT_LE(size, mat.capacity().first);
    EXPECT_LE(size, mat.capacity().second);
    EXPECT_NE(nullptr, mat.ptr());
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<std::complex<double>, dca::linalg::GPU> mat(name, size, capacity);
    EXPECT_EQ(name, mat.get_name());
    EXPECT_EQ(std::make_pair(size, size), mat.size());
    EXPECT_LE(capacity, mat.capacity().first);
    EXPECT_LE(capacity, mat.capacity().second);
    EXPECT_NE(nullptr, mat.ptr());
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));
  }
  {
    dca::linalg::Matrix<std::complex<float>, dca::linalg::GPU> mat(name, size2);
    EXPECT_EQ(name, mat.get_name());
    EXPECT_EQ(size2, mat.size());
    EXPECT_LE(size2.first, mat.capacity().first);
    EXPECT_LE(size2.second, mat.capacity().second);
    EXPECT_NE(nullptr, mat.ptr());
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));
  }
}

TEST(MatrixGPUTest, Properties) {
  {
    std::pair<int, int> size2(3, 5);
    std::pair<int, int> capacity2(5, 5);

    dca::linalg::Matrix<float, dca::linalg::GPU> mat(size2, capacity2);
    EXPECT_FALSE(mat.is_square());
    EXPECT_EQ(size2.first, mat.nrRows());
    EXPECT_EQ(size2.second, mat.nrCols());
    EXPECT_EQ(mat.capacity().first, mat.leadingDimension());
  }
  {
    std::pair<int, int> size2(5, 5);
    std::pair<int, int> capacity2(256, 5);

    dca::linalg::Matrix<float, dca::linalg::GPU> mat(size2, capacity2);
    EXPECT_TRUE(mat.is_square());
    EXPECT_EQ(size2.first, mat.nrRows());
    EXPECT_EQ(size2.second, mat.nrCols());
    EXPECT_EQ(mat.capacity().first, mat.leadingDimension());
  }
}

TEST(MatrixGPUTest, ElementPointers) {
  // Check if the pointers are computed correctly.
  std::pair<int, int> size2(5, 3);

  dca::linalg::Matrix<int, dca::linalg::GPU> mat(size2);
  const dca::linalg::Matrix<int, dca::linalg::GPU>& mat_const_ref(mat);
  for (int j = 0; j < mat.nrCols(); ++j)
    for (int i = 0; i < mat.nrRows(); ++i) {
      int* ptr = mat.ptr();
      int diff_ptr = i + j * mat.leadingDimension();
      EXPECT_EQ(diff_ptr, mat.ptr(i, j) - ptr);
      EXPECT_EQ(mat.ptr(i, j), mat_const_ref.ptr(i, j));
    }
}

TEST(MatrixGPUTest, CopyConstructor) {
  std::pair<int, int> size2(2, 3);

  dca::linalg::Matrix<float, dca::linalg::GPU> mat("name", size2);
  auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
  testing::setMatrixElements(mat, el_value);

  dca::linalg::Matrix<float, dca::linalg::GPU> mat_copy(mat, "another name");
  EXPECT_EQ("another name", mat_copy.get_name());
  EXPECT_EQ(mat.size(), mat_copy.size());
  EXPECT_LE(mat.size().first, mat_copy.capacity().first);
  EXPECT_LE(mat.size().second, mat_copy.capacity().second);

  for (int j = 0; j < mat.nrCols(); ++j)
    for (int i = 0; i < mat.nrRows(); ++i) {
      EXPECT_EQ(testing::getFromDevice(mat.ptr(i, j)), testing::getFromDevice(mat_copy.ptr(i, j)));
      EXPECT_NE(mat.ptr(i, j), mat_copy.ptr(i, j));
    }
}

TEST(MatrixGPUTest, MoveConstructor) {
  using MatrixType = dca::linalg::Matrix<double, dca::linalg::GPU>;
  MatrixType mat("matrix name", 4);
  auto el_value = [](int i, int j) { return 3.14 * i - 2.5 * j; };
  testing::setMatrixElements(mat, el_value);
  MatrixType mat_copy(mat);

  MatrixType thief(std::move(mat), "thief matrix");
  EXPECT_EQ("thief matrix", thief.get_name());

  for (int j = 0; j < mat.nrCols(); ++j)
    for (int i = 0; i < mat.nrRows(); ++i)
      EXPECT_EQ(testing::getFromDevice(mat_copy.ptr(i, j)), testing::getFromDevice(thief.ptr(i, j)));

  // The original matrix is empty.
  EXPECT_EQ(std::make_pair(0, 0), mat.size());
}

TEST(MatrixGPUTest, Assignement) {
  {
    // Assign a matrix that fits into the capacity.
    std::pair<int, int> size2(2, 3);

    dca::linalg::Matrix<float, dca::linalg::GPU> mat_copy(10);
    auto old_ptr = mat_copy.ptr();
    auto capacity = mat_copy.capacity();

    dca::linalg::Matrix<float, dca::linalg::GPU> mat("name", size2);
    auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    mat_copy = mat;
    EXPECT_EQ(mat.size(), mat_copy.size());
    EXPECT_EQ(capacity, mat_copy.capacity());
    EXPECT_EQ(old_ptr, mat_copy.ptr());

    for (int j = 0; j < mat.nrCols(); ++j)
      for (int i = 0; i < mat.nrRows(); ++i) {
        EXPECT_EQ(testing::getFromDevice(mat.ptr(i, j)), testing::getFromDevice(mat_copy.ptr(i, j)));
        EXPECT_NE(mat.ptr(i, j), mat_copy.ptr(i, j));
      }
  }
  {
    // Assign a matrix that does not fit into the capacity.
    dca::linalg::Matrix<float, dca::linalg::GPU> mat_copy(10);
    auto size2 = mat_copy.capacity();
    ++size2.first;

    dca::linalg::Matrix<float, dca::linalg::GPU> mat("name", size2);
    auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    mat_copy = mat;
    EXPECT_EQ(mat.size(), mat_copy.size());
    EXPECT_LE(mat.size().first, mat_copy.capacity().first);
    EXPECT_LE(mat.size().second, mat_copy.capacity().second);

    for (int j = 0; j < mat.nrCols(); ++j)
      for (int i = 0; i < mat.nrRows(); ++i) {
        EXPECT_EQ(testing::getFromDevice(mat.ptr(i, j)), testing::getFromDevice(mat_copy.ptr(i, j)));
        EXPECT_NE(mat.ptr(i, j), mat_copy.ptr(i, j));
      }
  }
}

TEST(MatrixGPUTest, MoveAssignement) {
  using MatrixType = dca::linalg::Matrix<float, dca::linalg::GPU>;
  MatrixType mat("matrix name", std::make_pair(2, 5));
  auto el_value = [](int i, int j) { return 3 * i + 2 * j; };
  testing::setMatrixElements(mat, el_value);
  MatrixType mat_copy(mat);

  MatrixType thief("thief name");
  thief = std::move(mat);

  EXPECT_EQ("thief name", thief.get_name());

  EXPECT_EQ(mat_copy.size(), thief.size());
  for (int j = 0; j < mat_copy.nrCols(); ++j)
    for (int i = 0; i < mat_copy.nrRows(); ++i)
      EXPECT_EQ(testing::getFromDevice(mat_copy.ptr(i, j)), testing::getFromDevice(thief.ptr(i, j)));

  // mat is now empty.
  EXPECT_EQ(std::make_pair(0, 0), mat.size());

  // Test return value.
  const MatrixType* const mat_ptr = &mat;
  EXPECT_EQ(mat_ptr, &(mat = std::move(thief)));
}

TEST(MatrixGPUTest, Set) {
  {
    // Assign a matrix that fits into the capacity.
    std::pair<int, int> size2(2, 3);

    dca::linalg::Matrix<float, dca::linalg::GPU> mat_copy(10);
    auto old_ptr = mat_copy.ptr();
    auto capacity = mat_copy.capacity();

    dca::linalg::Matrix<float, dca::linalg::GPU> mat("name", size2);
    auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    mat_copy.set(mat, 0, 1);
    EXPECT_EQ(mat.size(), mat_copy.size());
    EXPECT_EQ(capacity, mat_copy.capacity());
    EXPECT_EQ(old_ptr, mat_copy.ptr());

    for (int j = 0; j < mat.nrCols(); ++j)
      for (int i = 0; i < mat.nrRows(); ++i) {
        EXPECT_EQ(testing::getFromDevice(mat.ptr(i, j)), testing::getFromDevice(mat_copy.ptr(i, j)));
        EXPECT_NE(mat.ptr(i, j), mat_copy.ptr(i, j));
      }
  }
  {
    // Assign a matrix that does not fit into the capacity.
    dca::linalg::Matrix<float, dca::linalg::GPU> mat_copy(10);
    auto size2 = mat_copy.capacity();
    ++size2.first;

    dca::linalg::Matrix<float, dca::linalg::GPU> mat("name", size2);
    auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    mat_copy.set(mat, 0, 1);
    EXPECT_EQ(mat.size(), mat_copy.size());
    EXPECT_LE(mat.size().first, mat_copy.capacity().first);
    EXPECT_LE(mat.size().second, mat_copy.capacity().second);

    for (int j = 0; j < mat.nrCols(); ++j)
      for (int i = 0; i < mat.nrRows(); ++i) {
        EXPECT_EQ(testing::getFromDevice(mat.ptr(i, j)), testing::getFromDevice(mat_copy.ptr(i, j)));
        EXPECT_NE(mat.ptr(i, j), mat_copy.ptr(i, j));
      }
  }
}

TEST(MatrixGPUTest, Swap) {
  std::string mat1_name = "name 1";
  std::pair<int, int> mat1_size(7, 8);
  dca::linalg::Matrix<float, dca::linalg::GPU> mat1(mat1_name, mat1_size);
  auto mat1_capacity = mat1.capacity();
  auto mat1_ptr = mat1.ptr();

  std::string mat2_name = "name 2";
  std::pair<int, int> mat2_size(2, 128);
  dca::linalg::Matrix<float, dca::linalg::GPU> mat2(mat2_name, mat2_size);
  auto mat2_capacity = mat2.capacity();
  auto mat2_ptr = mat2.ptr();

  mat1.swap(mat2);
  EXPECT_EQ(mat1_name, mat1.get_name());
  EXPECT_EQ(mat1_size, mat2.size());
  EXPECT_EQ(mat1_capacity, mat2.capacity());
  EXPECT_EQ(mat1_ptr, mat2.ptr());

  EXPECT_EQ(mat2_name, mat2.get_name());
  EXPECT_EQ(mat2_size, mat1.size());
  EXPECT_EQ(mat2_capacity, mat1.capacity());
  EXPECT_EQ(mat2_ptr, mat1.ptr());
}

TEST(MatrixCPUTest, SwapWithName) {
  std::string mat1_name = "name 1";
  std::pair<int, int> mat1_size(7, 8);
  dca::linalg::Matrix<float, dca::linalg::GPU> mat1(mat1_name, mat1_size);
  auto mat1_capacity = mat1.capacity();
  auto mat1_ptr = mat1.ptr();

  std::string mat2_name = "name 2";
  std::pair<int, int> mat2_size(2, 128);
  dca::linalg::Matrix<float, dca::linalg::GPU> mat2(mat2_name, mat2_size);
  auto mat2_capacity = mat2.capacity();
  auto mat2_ptr = mat2.ptr();

  mat1.swapWithName(mat2);
  EXPECT_EQ(mat1_name, mat2.get_name());
  EXPECT_EQ(mat1_size, mat2.size());
  EXPECT_EQ(mat1_capacity, mat2.capacity());
  EXPECT_EQ(mat1_ptr, mat2.ptr());

  EXPECT_EQ(mat2_name, mat1.get_name());
  EXPECT_EQ(mat2_size, mat1.size());
  EXPECT_EQ(mat2_capacity, mat1.capacity());
  EXPECT_EQ(mat2_ptr, mat1.ptr());
}

TEST(MatrixGPUTest, ResizePair) {
  {
    std::pair<int, int> size2(4, 2);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);
    auto el_value = [](int i, int j) { return 1 + 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    // Resize to capacity. No reallocation has to take place.
    auto old_ptr = mat.ptr();
    auto capacity = mat.capacity();
    auto new_size = capacity;
    mat.resize(new_size);
    EXPECT_EQ(new_size, mat.size());
    EXPECT_EQ(capacity, mat.capacity());
    EXPECT_EQ(old_ptr, mat.ptr());

    // Check the value of the elements.
    for (int j = 0; j < size2.second; ++j)
      for (int i = 0; i < size2.first; ++i) {
        long el = el_value(i, j);
        EXPECT_EQ(el, testing::getFromDevice(mat.ptr(i, j)));
      }
  }
  {
    std::pair<int, int> size2(5, 2);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);
    auto old_ptr = mat.ptr();
    auto capacity = mat.capacity();
    auto el_value = [](int i, int j) { return 1 + 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    // Shrink the matrix. No reallocation has to take place.
    auto new_size = mat.size();
    --new_size.first;
    mat.resize(new_size);
    EXPECT_EQ(new_size, mat.size());
    EXPECT_EQ(capacity, mat.capacity());
    EXPECT_EQ(old_ptr, mat.ptr());

    // Check the value of the elements.
    for (int j = 0; j < mat.nrCols(); ++j)
      for (int i = 0; i < mat.nrRows(); ++i) {
        long el = el_value(i, j);
        EXPECT_EQ(el, testing::getFromDevice(mat.ptr(i, j)));
      }
  }
  {
    // We create a matrix with (padding_col + 1) colums,
    // where padding_col is the column capacity of a (1, 1) matrix.
    // In this way the column capacity decreases when we resize the matrix to 1 column.
    int padding_col =
        dca::linalg::Matrix<long, dca::linalg::GPU>(std::make_pair(1, 1)).capacity().second;
    std::pair<int, int> size2(5, padding_col + 1);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);
    auto old_ptr = mat.ptr();
    auto capacity = mat.capacity();
    auto el_value = [](int i, int j) { return 1 + 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    // The new number of rows is larger than capacity().first.
    // The number of columns is shrinking and the new column capacity will be smaller.
    // Reallocation has to take place.
    auto new_size = std::make_pair(capacity.first + 1, 1);
    mat.resize(new_size);
    EXPECT_EQ(new_size, mat.size());
    EXPECT_LE(new_size.first, mat.capacity().first);
    EXPECT_LE(new_size.second, mat.capacity().second);
    EXPECT_NE(old_ptr, mat.ptr());
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));

    // Check the value of the elements.
    for (int j = 0; j < 1; ++j)
      for (int i = 0; i < size2.first; ++i) {
        long el = el_value(i, j);
        EXPECT_EQ(el, testing::getFromDevice(mat.ptr(i, j)));
      }
  }
  {
    // We create a matrix with (padding_row + 1) colums,
    // where padding_row is the row capacity of a (1, 1) matrix.
    // In this way the row capacity decreases when we resize the matrix to 1 row.
    int padding_row =
        dca::linalg::Matrix<long, dca::linalg::GPU>(std::make_pair(1, 1)).capacity().first;
    std::pair<int, int> size2(padding_row + 1, 2);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);
    auto old_ptr = mat.ptr();
    auto capacity = mat.capacity();
    auto el_value = [](int i, int j) { return 1 + 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    // The new number of columns is larger than capacity().second.
    // The number of rows is shrinking and the new row capacity will be smaller.
    // Reallocation has to take place.
    auto new_size = std::make_pair(1, capacity.second + 1);
    mat.resize(new_size);
    EXPECT_EQ(new_size, mat.size());
    EXPECT_LE(new_size.first, mat.capacity().first);
    EXPECT_LE(new_size.second, mat.capacity().second);
    EXPECT_NE(old_ptr, mat.ptr());
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));

    // Check the value of the elements.
    for (int j = 0; j < size2.second; ++j)
      for (int i = 0; i < 1; ++i) {
        long el = el_value(i, j);
        EXPECT_EQ(el, testing::getFromDevice(mat.ptr(i, j)));
      }
  }
}

TEST(MatrixGPUTest, ResizeValue) {
  {
    std::pair<int, int> size2(4, 2);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);

    auto el_value = [](int i, int j) { return 1 + 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    // Resize to capacity. No reallocation has to take place.
    auto old_ptr = mat.ptr();
    auto capacity = mat.capacity();
    int new_size = std::min(capacity.first, capacity.second);
    mat.resize(new_size);
    EXPECT_EQ(std::make_pair(new_size, new_size), mat.size());
    EXPECT_EQ(capacity, mat.capacity());
    EXPECT_EQ(old_ptr, mat.ptr());

    // Check the value of the elements.
    for (int j = 0; j < mat.nrCols(); ++j)
      for (int i = 0; i < mat.nrRows(); ++i) {
        long el = el_value(i, j);
        EXPECT_EQ(el, testing::getFromDevice(mat.ptr(i, j)));
      }
  }
  {
    std::pair<int, int> size2(5, 3);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);
    auto old_ptr = mat.ptr();
    auto capacity = mat.capacity();
    auto el_value = [](int i, int j) { return 1 + 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    // Shrink the matrix. No reallocation has to take place.
    int new_size = 2;
    mat.resize(new_size);
    EXPECT_EQ(std::make_pair(new_size, new_size), mat.size());
    EXPECT_EQ(capacity, mat.capacity());
    EXPECT_EQ(old_ptr, mat.ptr());

    // Check the value of the elements.
    for (int j = 0; j < mat.nrCols(); ++j)
      for (int i = 0; i < mat.nrRows(); ++i) {
        long el = el_value(i, j);
        EXPECT_EQ(el, testing::getFromDevice(mat.ptr(i, j)));
      }
  }
  {
    std::pair<int, int> size2(3, 3);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);
    auto old_ptr = mat.ptr();
    auto capacity = mat.capacity();
    auto el_value = [](int i, int j) { return 1 + 3 * i - 2 * j; };
    testing::setMatrixElements(mat, el_value);

    // New size is larger than capacity.
    // Reallocation has to take place.
    int new_size = std::min(capacity.first, capacity.second) + 1;
    mat.resize(new_size);
    EXPECT_EQ(std::make_pair(new_size, new_size), mat.size());
    EXPECT_LE(new_size, mat.capacity().first);
    EXPECT_LE(new_size, mat.capacity().second);
    EXPECT_NE(old_ptr, mat.ptr());
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));

    // Check the value of the elements.
    for (int j = 0; j < size2.second; ++j)
      for (int i = 0; i < size2.first; ++i) {
        long el = el_value(i, j);
        EXPECT_EQ(el, testing::getFromDevice(mat.ptr(i, j)));
      }
  }
}

TEST(MatrixGPUTest, ResizeNoCopyPair) {
  {
    std::pair<int, int> size2(4, 2);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);

    // Resize to capacity. No reallocation has to take place.
    auto old_ptr = mat.ptr();
    auto capacity = mat.capacity();
    auto new_size = capacity;
    mat.resizeNoCopy(new_size);
    EXPECT_EQ(new_size, mat.size());
    EXPECT_EQ(capacity, mat.capacity());
    EXPECT_EQ(old_ptr, mat.ptr());
  }
  {
    std::pair<int, int> size2(5, 2);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);
    auto old_ptr = mat.ptr();
    auto capacity = mat.capacity();

    // Shrink the matrix. No reallocation has to take place.
    auto new_size = mat.size();
    --new_size.first;
    mat.resizeNoCopy(new_size);
    EXPECT_EQ(new_size, mat.size());
    EXPECT_EQ(capacity, mat.capacity());
    EXPECT_EQ(old_ptr, mat.ptr());
  }
  {
    std::pair<int, int> size2(5, 2);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);
    auto capacity = mat.capacity();

    // New number of rows is larger than capacity().first.
    // Reallocation has to take place.
    auto new_size = std::make_pair(capacity.first + 1, 1);
    mat.resizeNoCopy(new_size);
    EXPECT_EQ(new_size, mat.size());
    EXPECT_LE(new_size.first, mat.capacity().first);
    EXPECT_LE(new_size.second, mat.capacity().second);
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));
  }
  {
    std::pair<int, int> size2(5, 2);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);
    auto capacity = mat.capacity();

    // New number of columns is larger than capacity().second.
    // Reallocation has to take place.
    auto new_size = std::make_pair(1, capacity.second + 1);
    mat.resizeNoCopy(new_size);
    EXPECT_EQ(new_size, mat.size());
    EXPECT_LE(new_size.first, mat.capacity().first);
    EXPECT_LE(new_size.second, mat.capacity().second);
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));
  }
}

TEST(MatrixGPUTest, ResizeNoCopyValue) {
  {
    std::pair<int, int> size2(4, 2);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);

    // Resize to capacity. No reallocation has to take place.
    auto old_ptr = mat.ptr();
    auto capacity = mat.capacity();
    int new_size = std::min(capacity.first, capacity.second);
    mat.resizeNoCopy(new_size);
    EXPECT_EQ(std::make_pair(new_size, new_size), mat.size());
    EXPECT_EQ(capacity, mat.capacity());
    EXPECT_EQ(old_ptr, mat.ptr());
  }
  {
    std::pair<int, int> size2(5, 3);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);
    auto old_ptr = mat.ptr();
    auto capacity = mat.capacity();

    // Shrink the matrix. No reallocation has to take place.
    int new_size = 2;
    mat.resizeNoCopy(new_size);
    EXPECT_EQ(std::make_pair(new_size, new_size), mat.size());
    EXPECT_EQ(capacity, mat.capacity());
    EXPECT_EQ(old_ptr, mat.ptr());
  }
  {
    std::pair<int, int> size2(3, 3);

    dca::linalg::Matrix<long, dca::linalg::GPU> mat(size2);
    auto capacity = mat.capacity();

    // New size is larger than capacity.
    // Reallocation has to take place.
    int new_size = std::min(capacity.first, capacity.second) + 1;
    mat.resizeNoCopy(new_size);
    EXPECT_EQ(std::make_pair(new_size, new_size), mat.size());
    EXPECT_LE(new_size, mat.capacity().first);
    EXPECT_LE(new_size, mat.capacity().second);
    EXPECT_TRUE(testing::isDevicePointer(mat.ptr()));
  }
}

TEST(MatrixGPUTest, Clear) {
  dca::linalg::Matrix<float, dca::linalg::GPU> mat(42);

  EXPECT_EQ(std::make_pair(42, 42), mat.size());
  mat.clear();
  EXPECT_EQ(std::make_pair(0, 0), mat.size());
  EXPECT_EQ(std::make_pair(0, 0), mat.capacity());
}

TEST(MatrixGPUTest, setToZero) {
  dca::linalg::Matrix<long, dca::linalg::GPU> mat(3);
  auto func = [](int i, int j) { return 10 * i - j; };
  testing::setMatrixElements(mat, func);

  cudaStream_t stream;
  cudaStreamCreate(&stream);
  mat.setToZero(stream);
  cudaStreamSynchronize(stream);

  dca::linalg::Matrix<long, dca::linalg::CPU> mat_copy(mat);

  for (int j = 0; j < mat_copy.nrCols(); ++j)
    for (int i = 0; i < mat_copy.nrRows(); ++i)
      EXPECT_EQ(0, mat_copy(i, j));

  cudaStreamDestroy(stream);
}
