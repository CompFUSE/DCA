// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the MatrixView<CPU> class.

#include "dca/linalg/matrix_view.hpp"

#include <complex>

#include "gtest/gtest.h"

#include "cpu_test_util.hpp"
#include "dca/linalg/matrix.hpp"

TEST(MatrixViewTest, Constructors) {
  std::pair<int, int> size(4, 5);
  {
    dca::linalg::Matrix<float, dca::linalg::CPU> mat(size);
    dca::linalg::MatrixView<float, dca::linalg::CPU> view(mat);
    EXPECT_EQ(mat.nrRows(), view.nrRows());
    EXPECT_EQ(mat.nrCols(), view.nrCols());
    EXPECT_EQ(mat.ptr(), view.ptr());
    EXPECT_EQ(mat.leadingDimension(), view.leadingDimension());
  }
  {
    dca::linalg::Matrix<double, dca::linalg::CPU> mat(size);
    const int delta_i(1), delta_j(2);
    dca::linalg::MatrixView<double, dca::linalg::CPU> view(mat, delta_i, delta_j);
    EXPECT_EQ(mat.nrRows() - delta_i, view.nrRows());
    EXPECT_EQ(mat.nrCols() - delta_j, view.nrCols());
    EXPECT_EQ(mat.ptr(delta_i, delta_j), view.ptr());
    EXPECT_EQ(mat.leadingDimension(), view.leadingDimension());
  }
  {
    dca::linalg::Matrix<double, dca::linalg::CPU> mat(size);
    const int delta_i(0), delta_j(3);
    const int ni(1), nj(0);
    dca::linalg::MatrixView<double, dca::linalg::CPU> view(mat, delta_i, delta_j, ni, nj);
    EXPECT_EQ(ni, view.nrRows());
    EXPECT_EQ(nj, view.nrCols());
    EXPECT_EQ(mat.ptr(delta_i, delta_j), view.ptr());
    EXPECT_EQ(mat.leadingDimension(), view.leadingDimension());
  }
}

TEST(MatrixViewTest, ReadWrite) {
  dca::linalg::Matrix<ushort, dca::linalg::CPU> mat(4);
  dca::linalg::MatrixView<ushort, dca::linalg::CPU> view(mat);

  view(1, 2) = 2;
  EXPECT_EQ(2, mat(1, 2));

  mat(2, 3) = 1;
  EXPECT_EQ(1, view(2, 3));

  dca::linalg::MatrixView<ushort, dca::linalg::CPU> view_shifted(mat, 1, 2);
  EXPECT_EQ(2, view_shifted(0, 0));

  EXPECT_DEBUG_DEATH(view(-1, 2), "Assertion.");
  EXPECT_DEBUG_DEATH(view(0, 4), "Assertion.");
}

TEST(MatrixViewTest, Assignment) {
  dca::linalg::Matrix<int, dca::linalg::CPU> mat(4);
  auto init_func = [](int i, int j) { return i >= 2 ? 1 : 0; };
  testing::setMatrixElements(mat, init_func);

  dca::linalg::MatrixView<int, dca::linalg::CPU> upper_left(mat, 0, 0, 2, 2);
  dca::linalg::MatrixView<int, dca::linalg::CPU> lower_right(mat, 2, 2, 2, 2);

  upper_left = lower_right;  // Assign ones to upper left submatrix.
  for (int j = 0; j < 2; ++j)
    for (int i = 0; i < 2; ++i)
      EXPECT_EQ(1, mat(i, j));

  dca::linalg::MatrixView<int, dca::linalg::CPU> another_size(mat, 0, 0, 3, 3);
  // Invalid assignment:
  EXPECT_THROW(upper_left = another_size       , std::invalid_argument);
}

TEST(MatrixViewTest, MakeConstantView) {
  dca::linalg::Matrix<double, dca::linalg::CPU> mat(std::make_pair(4, 2));
  auto init_func = [](int i, int j) { return j + 10. * i; };
  testing::setMatrixElements(mat, init_func);
  {
    const dca::linalg::Matrix<double, dca::linalg::CPU> const_mat(mat);
    auto const_view_ptr = dca::linalg::makeConstantView(const_mat);
    const auto& const_view = *const_view_ptr;
    EXPECT_EQ(const_mat.nrRows(), const_view.nrRows());
    EXPECT_EQ(const_mat.nrCols(), const_view.nrCols());
    EXPECT_EQ(const_mat.ptr(), const_view.ptr());
    EXPECT_EQ(const_mat.leadingDimension(), const_view.leadingDimension());

    for (int j = 0; j < const_view.nrCols(); ++j)
      for (int i = 0; i < const_view.nrCols(); ++i)
        EXPECT_EQ(const_mat(i, j), const_view(i, j));
  }
  {
    const dca::linalg::Matrix<double, dca::linalg::CPU> const_mat(mat);
    auto const_view_ptr = dca::linalg::makeConstantView(const_mat, 1, 0);
    const auto& const_view = *const_view_ptr;
    EXPECT_EQ(const_mat.nrRows() - 1, const_view.nrRows());
    EXPECT_EQ(const_mat.nrCols(), const_view.nrCols());
    EXPECT_EQ(const_mat.ptr(1, 0), const_view.ptr());
    EXPECT_EQ(const_mat.leadingDimension(), const_view.leadingDimension());

    for (int j = 0; j < const_view.nrCols(); ++j)
      for (int i = 0; i < const_view.nrCols(); ++i)
        EXPECT_EQ(const_mat(i + 1, j), const_view(i, j));
  }
  {
    const dca::linalg::Matrix<double, dca::linalg::CPU> const_mat(mat);
    auto const_view_ptr = dca::linalg::makeConstantView(const_mat, 0, 1, 2, 1);
    const auto& const_view = *const_view_ptr;
    EXPECT_EQ(2, const_view.nrRows());
    EXPECT_EQ(1, const_view.nrCols());
    EXPECT_EQ(const_mat.ptr(0, 1), const_view.ptr());
    EXPECT_EQ(const_mat.leadingDimension(), const_view.leadingDimension());

    for (int j = 0; j < const_view.nrCols(); ++j)
      for (int i = 0; i < const_view.nrCols(); ++i)
        EXPECT_EQ(const_mat(i, j + 1), const_view(i, j));
  }
}
