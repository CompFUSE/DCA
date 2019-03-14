// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the ReshapableMatrix<GPU> class.

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
    dca::linalg::ReshapableMatrix<int, dca::linalg::GPU> mat;
    EXPECT_LE(0, mat.nrRows());
    EXPECT_LE(0, mat.nrCols());
    EXPECT_EQ(0, mat.leadingDimension());
    EXPECT_LE(0, mat.capacity());
    EXPECT_EQ(nullptr, mat.ptr());
  }
  {
    const int size = 5;
    dca::linalg::ReshapableMatrix<double, dca::linalg::GPU> mat(size);
    ASSERT_EQ(std::make_pair(size, size), mat.size());
    ASSERT_EQ(size, mat.nrRows());
    ASSERT_EQ(size, mat.nrCols());
    EXPECT_EQ(size, mat.leadingDimension());
    ASSERT_LE(size, mat.capacity());
    ASSERT_NE(nullptr, mat.ptr());
  }
  {
    const auto size = std::make_pair(3, 2);
    dca::linalg::ReshapableMatrix<float, dca::linalg::GPU> mat(size);
    ASSERT_EQ(size, mat.size());
    ASSERT_EQ(size.first, mat.nrRows());
    ASSERT_EQ(size.second, mat.nrCols());
    EXPECT_EQ(size.first, mat.leadingDimension());
    ASSERT_LE(product(size), mat.capacity());
    ASSERT_NE(nullptr, mat.ptr());
  }
}

TEST(MatrixCPUTest, CopyConstructor) {
  dca::linalg::ReshapableMatrix<float, dca::linalg::CPU> mat_cpu(std::make_pair(5, 8));
  auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
  testing::setMatrixElements(mat_cpu, el_value);

  dca::linalg::ReshapableMatrix<float, dca::linalg::GPU> mat_gpu(mat_cpu);
  EXPECT_EQ(mat_gpu.size(), mat_cpu.size());

  dca::linalg::ReshapableMatrix<float, dca::linalg::GPU> mat_gpu2(mat_gpu);
  EXPECT_EQ(mat_gpu2.size(), mat_gpu.size());

  dca::linalg::ReshapableMatrix<float, dca::linalg::CPU> mat_cpu2(mat_gpu2);

  EXPECT_EQ(mat_cpu, mat_cpu2);
}

TEST(MatrixCPUTest, MoveConstructor) {
  dca::linalg::ReshapableMatrix<float, dca::linalg::GPU> mat(std::make_pair(5, 8));
  const float* mat_ptr = mat.ptr();
  auto size = mat.size();
  std::size_t capacity = mat.capacity();

  dca::linalg::ReshapableMatrix<float, dca::linalg::GPU> mat2(std::move(mat));

  EXPECT_EQ(size, mat2.size());
  EXPECT_EQ(capacity, mat2.capacity());
  EXPECT_EQ(mat_ptr, mat2.ptr());
}

TEST(MatrixCPUTest, Assignement) {
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> mat_cpu(std::make_pair(3, 15));
  dca::linalg::Matrix<std::complex<double>, dca::linalg::GPU> mat_gpu(10);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::GPU> mat_gpu2;

  auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
  testing::setMatrixElements(mat_cpu, el_value);

  // Test chain assignment.
  mat_gpu2 = mat_gpu = mat_cpu;

  EXPECT_EQ(mat_cpu.size(), mat_gpu.size());
  EXPECT_EQ(mat_cpu.size(), mat_gpu2.size());

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> mat_cpu2(mat_gpu2);
  EXPECT_EQ(mat_cpu2, mat_cpu);
}

TEST(MatrixCPUTest, MoveAssignement) {
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> mat_cpu(6);

  auto el_value = [](int i, int j) { return 3 * i - 2 * j; };
  testing::setMatrixElements(mat_cpu, el_value);

  dca::linalg::Matrix<std::complex<double>, dca::linalg::GPU> mat_gpu(mat_cpu);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::GPU> mat_gpu2;

  mat_gpu2 = std::move(mat_gpu);

  EXPECT_EQ(mat_cpu.size(), mat_gpu2.size());

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> mat_cpu2(mat_gpu2);

  EXPECT_EQ(mat_cpu2, mat_cpu);
}

TEST(MatrixCPUTest, SetAsync) {
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> mat_cpu(8);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> mat_cpu_out(8);
  dca::linalg::Matrix<std::complex<double>, dca::linalg::GPU> mat_gpu(8);

  testing::setMatrixElements(mat_cpu, [](int i, int j) { return 3 * i * i - j; });

  dca::linalg::util::CudaStream stream;

  mat_gpu.setAsync(mat_cpu, stream);
  mat_cpu_out.setAsync(mat_gpu, stream);
  cudaStreamSynchronize(stream);

  EXPECT_EQ(mat_cpu, mat_cpu_out);
}

TEST(MatrixCPUTest, ResizeNoCopy) {
  dca::linalg::ReshapableMatrix<int, dca::linalg::GPU, dca::linalg::util::ManagedAllocator<int>> mat(
      std::make_pair(3, 7));

  const int* old_ptr = mat.ptr();

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

TEST(MatrixCPUTest, Clear) {
  dca::linalg::ReshapableMatrix<float, dca::linalg::GPU> mat(42);

  EXPECT_EQ(std::make_pair(42, 42), mat.size());
  mat.clear();
  EXPECT_EQ(std::make_pair(0, 0), mat.size());
  EXPECT_EQ(0, mat.capacity());
}
