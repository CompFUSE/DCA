// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the ReshapableMatrix<CPU> class.

#include "dca/linalg/reshapable_matrix.hpp"
#include <complex>
#include <string>
#include "gtest/gtest.h"
#include "cpu_test_util.hpp"

TEST(MatrixCPUTest, Constructor) {
  std::pair<int, int> size_pair(4, 5);
  std::pair<int, int> capacity2(13, 17);

  auto product = [](auto& p) { return p.first * p.second; };

  {
    dca::linalg::ReshapableMatrix<int, dca::linalg::CPU> mat;
    EXPECT_LE(0, mat.nrRows());
    EXPECT_LE(0, mat.nrCols());
    EXPECT_EQ(0, mat.leadingDimension());
    EXPECT_LE(0, mat.capacity());
    EXPECT_EQ(nullptr, mat.ptr());
  }
  {
    const int size = 5;
    dca::linalg::ReshapableMatrix<double, dca::linalg::CPU> mat(size);
    ASSERT_EQ(std::make_pair(size, size), mat.size());
    ASSERT_EQ(size, mat.nrRows());
    ASSERT_EQ(size, mat.nrCols());
    EXPECT_EQ(size, mat.leadingDimension());
    ASSERT_LE(size, mat.capacity());
    ASSERT_NE(nullptr, mat.ptr());
  }
  {
    const auto size = std::make_pair(3, 2);
    dca::linalg::ReshapableMatrix<float, dca::linalg::CPU> mat(size);
    ASSERT_EQ(size, mat.size());
    ASSERT_EQ(size.first, mat.nrRows());
    ASSERT_EQ(size.second, mat.nrCols());
    EXPECT_EQ(size.first, mat.leadingDimension());
    ASSERT_LE(product(size), mat.capacity());
    ASSERT_NE(nullptr, mat.ptr());
  }
}

TEST(MatrixCPUTest, CopyConstructor) {
  std::pair<int, int> size2(2, 3);
  dca::linalg::ReshapableMatrix<float, dca::linalg::CPU> mat(size2);
  auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
  testing::setMatrixElements(mat, el_value);

  dca::linalg::ReshapableMatrix<float, dca::linalg::CPU> mat_copy(mat);
  EXPECT_EQ(mat, mat_copy);
  EXPECT_LE(mat.size().first * mat.size().second, mat_copy.capacity());

  EXPECT_NE(mat.ptr(), mat_copy.ptr());
}

TEST(MatrixCPUTest, MoveConstructor) {
  using MatrixType = dca::linalg::ReshapableMatrix<double, dca::linalg::CPU>;
  std::pair<int, int> size(2, 3);
  MatrixType mat(size);

  auto f = [](int i, int j) { return 3.14 * i - 2.5 * j; };
  testing::setMatrixElements(mat, f);
  MatrixType mat_copy(mat);

  MatrixType mat_thief(std::move(mat));
  EXPECT_EQ(mat_copy, mat_thief);
  // The original matrix is now empty.
  EXPECT_EQ(std::make_pair(0, 0), mat.size());
}

TEST(MatrixCPUTest, Assignement) {
  dca::linalg::Matrix<float, dca::linalg::CPU> mat(6);
  dca::linalg::Matrix<float, dca::linalg::CPU> mat2(10);

  auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
  testing::setMatrixElements(mat2, el_value);

  mat = mat2;
  EXPECT_EQ(mat, mat2);
  EXPECT_NE(mat.ptr(), mat2.ptr());
}

TEST(MatrixCPUTest, MoveAssignement) {
  dca::linalg::Matrix<float, dca::linalg::CPU> mat(6);
  dca::linalg::Matrix<float, dca::linalg::CPU> mat2(10);

  auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
  testing::setMatrixElements(mat2, el_value);
  dca::linalg::Matrix<float, dca::linalg::CPU> mat2_copy(mat2);

  const float* mat2_ptr = mat2.ptr();

  mat = std::move(mat2);
  EXPECT_EQ(mat, mat2_copy);
  EXPECT_EQ(mat.ptr(), mat2_ptr);

  // mat2 was moved.
  EXPECT_NE(mat, mat2);
}

TEST(MatrixCPUTest, ComparisonOperators) {
  dca::linalg::ReshapableMatrix<float, dca::linalg::CPU> mat1(std::make_pair(2, 2));
  dca::linalg::ReshapableMatrix<float, dca::linalg::CPU> mat2(std::make_pair(2, 2));
  dca::linalg::ReshapableMatrix<float, dca::linalg::CPU> mat3(std::make_pair(2, 3));

  auto fill_func = [](int i, int j) { return std::cos(i) + j; };
  testing::setMatrixElements(mat1, fill_func);
  testing::setMatrixElements(mat2, fill_func);
  testing::setMatrixElements(mat3, fill_func);
  mat2(0, 1) = -1;

  // Matrices with different elements are not equal.
  EXPECT_TRUE(mat1 != mat2);
  EXPECT_FALSE(mat1 == mat2);
  // Matrices with different shapes are not equal.
  EXPECT_TRUE(mat1 != mat3);
  EXPECT_FALSE(mat1 == mat3);

  // Matrices with no elements are considered equal.
  mat1.resizeNoCopy(std::make_pair(0, 2));
  mat3.resizeNoCopy(std::make_pair(3, 0));
  EXPECT_TRUE(mat1 == mat3);
  EXPECT_FALSE(mat1 != mat3);
}

TEST(MatrixCPUTest, ElementPointers) {
  // Check if the pointers are computed correctly.
  std::pair<int, int> size2(5, 3);

  dca::linalg::ReshapableMatrix<int, dca::linalg::CPU> mat(size2);
  const dca::linalg::ReshapableMatrix<int, dca::linalg::CPU>& mat_const_ref(mat);
  for (int j = 0; j < mat.nrCols(); ++j)
    for (int i = 0; i < mat.nrRows(); ++i) {
      int* ptr = mat.ptr();
      int diff_ptr = i + j * mat.leadingDimension();
      EXPECT_EQ(diff_ptr, mat.ptr(i, j) - ptr);
      EXPECT_EQ(mat.ptr(i, j), mat_const_ref.ptr(i, j));
      EXPECT_EQ(mat.ptr(i, j), &mat(i, j));
      EXPECT_EQ(mat.ptr(i, j), &mat_const_ref(i, j));
    }
}

TEST(MatrixCPUTest, ElementAccess) {
  // Check if the different element accesses return the same value.
  std::pair<int, int> size2(2, 3);

  dca::linalg::ReshapableMatrix<int, dca::linalg::CPU> mat(size2);
  const dca::linalg::ReshapableMatrix<int, dca::linalg::CPU>& mat_const_ref(mat);
  for (int j = 0; j < mat.nrCols(); ++j)
    for (int i = 0; i < mat.nrRows(); ++i) {
      int el = 3 * i - 2 * j;
      mat(i, j) = el;
      EXPECT_EQ(el, mat(i, j));
      EXPECT_EQ(el, mat_const_ref(i, j));
      EXPECT_EQ(el, *(mat.ptr(i, j)));
      EXPECT_EQ(el, *(mat_const_ref.ptr(i, j)));
    }
}

TEST(MatrixCPUTest, ResizeNoCopy) {
  dca::linalg::ReshapableMatrix<std::complex<float>, dca::linalg::CPU> mat(std::make_pair(3, 7));
  const std::complex<float>* old_ptr = mat.ptr();

  auto product = [](std::pair<int, int>& p) { return p.first * p.second; };

  {
    // Resize inside the capacity.
    auto new_size = std::make_pair(7, 2);
    EXPECT_GE(mat.capacity(), product(new_size));

    const bool reallocation = mat.resizeNoCopy(new_size);
    EXPECT_FALSE(reallocation);
    EXPECT_EQ(new_size, mat.size());
    EXPECT_EQ(7, mat.leadingDimension());
    EXPECT_EQ(old_ptr, mat.ptr());
  }
  {
    // Resize outside the capacity.
    auto new_size = std::make_pair(17, 5);
    EXPECT_GE(product(new_size), mat.capacity());

    const bool reallocation = mat.resizeNoCopy(new_size);
    EXPECT_TRUE(reallocation);
    EXPECT_EQ(new_size, mat.size());
    EXPECT_EQ(17, mat.leadingDimension());
  }
}

TEST(MatrixCPUTest, SetAsync) {
  dca::linalg::ReshapableMatrix<float, dca::linalg::CPU> mat(8);
  dca::linalg::ReshapableMatrix<float, dca::linalg::CPU> mat2(8);

  testing::setMatrixElements(mat, [](int i, int j) { return 3 * i * i - j; });

  mat2.setAsync(mat, 0, 0);

  EXPECT_EQ(mat, mat2);
}

TEST(MatrixCPUTest, Clear) {
  dca::linalg::ReshapableMatrix<float, dca::linalg::CPU> mat(42);

  EXPECT_EQ(std::make_pair(42, 42), mat.size());
  mat.clear();
  EXPECT_EQ(std::make_pair(0, 0), mat.size());
  EXPECT_EQ(0, mat.capacity());
}
