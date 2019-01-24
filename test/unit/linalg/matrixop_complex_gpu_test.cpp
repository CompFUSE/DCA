// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the matrix operations with the Matrix<GPU> class with complex element type.

#include "dca/linalg/matrixop.hpp"
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include "gtest/gtest.h"
#include "dca/linalg/blas/blas3.hpp"
#include "dca/linalg/lapack/use_device.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/util/util_cublas.hpp"
#include "cpu_test_util.hpp"
#include "gpu_test_util.hpp"

template <typename ComplexType>
class MatrixopComplexGPUTest : public ::testing::Test {
public:
  static const typename ComplexType::value_type epsilon;
};
template <typename ComplexType>
const typename ComplexType::value_type MatrixopComplexGPUTest<ComplexType>::epsilon =
    std::numeric_limits<typename ComplexType::value_type>::epsilon();

typedef ::testing::Types<std::complex<float>, std::complex<double>> ComplexFloatingPointTypes;
TYPED_TEST_CASE(MatrixopComplexGPUTest, ComplexFloatingPointTypes);

TYPED_TEST(MatrixopComplexGPUTest, Gemm) {
  using ScalarType = TypeParam;
  auto val_a = [](int i, int j) { return ScalarType(3 * i - 2 * j, 1 - i * i + j); };
  auto val_b = [](int i, int j) { return ScalarType(i - 2 * j * j, 4 * i - 3 * j); };
  auto val_c = [](int i, int j) { return ScalarType(2 * i - j, 2 * j - i); };
  {
    std::pair<int, int> size_a(2, 3);
    std::pair<int, int> size_b(3, 4);
    std::pair<int, int> size_c(2, 4);

    {
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size_b);
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(size_c);

      testing::setMatrixElements(a, val_a);
      testing::setMatrixElements(b, val_b);
      testing::setMatrixElements(c, val_c);

      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> da(a);
      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> db(b);
      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(c);

      // Use CPU version as reference
      dca::linalg::matrixop::gemm(a, b, c);
      dca::linalg::matrixop::gemm(da, db, dc);

      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c_test(dc);

      for (int j = 0; j < c.nrCols(); ++j)
        for (int i = 0; i < c.nrRows(); ++i) {
          EXPECT_GE(500 * this->epsilon, std::abs(c(i, j) - c_test(i, j)));
        }
    }
    {
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size_b);
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(size_c);

      testing::setMatrixElements(a, val_a);
      testing::setMatrixElements(b, val_b);
      testing::setMatrixElements(c, val_c);

      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> da(a);
      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> db(b);
      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(c);

      ScalarType alpha = .57;
      ScalarType beta = -.712;

      // Use CPU version as reference
      dca::linalg::matrixop::gemm(alpha, a, b, beta, c);
      dca::linalg::matrixop::gemm(alpha, da, db, beta, dc);

      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c_test(dc);

      for (int j = 0; j < c.nrCols(); ++j)
        for (int i = 0; i < c.nrRows(); ++i) {
          EXPECT_GE(500 * this->epsilon, std::abs(c(i, j) - c_test(i, j)));
        }
    }
  }
  {
    int m = 2;
    int k = 4;
    int n = 3;
    std::pair<int, int> size_c(2, 3);

    char trans_opts[] = {'N', 'T', 'C'};

    for (auto transa : trans_opts) {
      std::pair<int, int> size_a;
      if (transa == 'N') {
        size_a.first = m;
        size_a.second = k;
      }
      else {
        size_a.first = k;
        size_a.second = m;
      }

      for (auto transb : trans_opts) {
        std::pair<int, int> size_b;
        if (transb == 'N') {
          size_b.first = k;
          size_b.second = n;
        }
        else {
          size_b.first = n;
          size_b.second = k;
        }

        {
          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size_b);
          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(size_c);

          testing::setMatrixElements(a, val_a);
          testing::setMatrixElements(b, val_b);
          testing::setMatrixElements(c, val_c);

          dca::linalg::Matrix<ScalarType, dca::linalg::GPU> da(a);
          dca::linalg::Matrix<ScalarType, dca::linalg::GPU> db(b);
          dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(c);

          // Use CPU version as reference
          dca::linalg::matrixop::gemm(transa, transb, a, b, c);
          dca::linalg::matrixop::gemm(transa, transb, da, db, dc);

          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c_test(dc);

          for (int j = 0; j < c.nrCols(); ++j)
            for (int i = 0; i < c.nrRows(); ++i) {
              EXPECT_GE(500 * this->epsilon, std::abs(c(i, j) - c_test(i, j)));
            }
        }
        {
          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size_b);
          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(size_c);

          testing::setMatrixElements(a, val_a);
          testing::setMatrixElements(b, val_b);
          testing::setMatrixElements(c, val_c);

          dca::linalg::Matrix<ScalarType, dca::linalg::GPU> da(a);
          dca::linalg::Matrix<ScalarType, dca::linalg::GPU> db(b);
          dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(c);

          ScalarType alpha = .57;
          ScalarType beta = -.712;

          // Use CPU version as reference
          dca::linalg::matrixop::gemm(transa, transb, alpha, a, b, beta, c);
          dca::linalg::matrixop::gemm(transa, transb, alpha, da, db, beta, dc);

          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c_test(dc);

          for (int j = 0; j < c.nrCols(); ++j)
            for (int i = 0; i < c.nrRows(); ++i) {
              EXPECT_GE(500 * this->epsilon, std::abs(c(i, j) - c_test(i, j)));
            }
        }
      }
    }
  }
}

TYPED_TEST(MatrixopComplexGPUTest, MultiplyDiagonal) {
  using ScalarType = TypeParam;
  std::pair<int, int> size_a(37, 45);
  auto val_a = [](int i, int j) { return ScalarType(3 * i - 2 * j, 1 - i * i + j); };
  auto val_d = [](int i) { return ScalarType(i, 1 - i); };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
  testing::setMatrixElements(a, val_a);
  dca::linalg::Matrix<ScalarType, dca::linalg::GPU> da(a);
  {
    dca::linalg::Vector<ScalarType, dca::linalg::CPU> d(a.nrRows());
    testing::setVectorElements(d, val_d);
    {
      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> db(size_a);

      // Test CPU vector.
      dca::linalg::matrixop::multiplyDiagonalLeft(d, da, db);
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(db);

      for (int j = 0; j < a.nrCols(); ++j)
        for (int i = 0; i < a.nrRows(); ++i) {
          EXPECT_GE(10 * this->epsilon, std::abs(b(i, j) - d[i] * a(i, j)));
        }
    }
    {
      dca::linalg::Vector<ScalarType, dca::linalg::GPU> dd(d);
      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> db(size_a);
      // Test GPU vector.

      dca::linalg::matrixop::multiplyDiagonalLeft(dd, da, db);
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(db);

      for (int j = 0; j < a.nrCols(); ++j)
        for (int i = 0; i < a.nrRows(); ++i) {
          EXPECT_GE(10 * this->epsilon, std::abs(b(i, j) - d[i] * a(i, j)));
        }
    }
  }
  {
    dca::linalg::Vector<ScalarType, dca::linalg::CPU> d(a.nrCols());
    testing::setVectorElements(d, val_d);
    {
      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> db(size_a);

      // Test CPU vector.
      dca::linalg::matrixop::multiplyDiagonalRight(da, d, db);
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(db);

      for (int j = 0; j < a.nrCols(); ++j)
        for (int i = 0; i < a.nrRows(); ++i) {
          EXPECT_GE(10 * this->epsilon, std::abs(b(i, j) - d[j] * a(i, j)));
        }
    }
    {
      dca::linalg::Vector<ScalarType, dca::linalg::GPU> dd(d);
      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> db(size_a);
      // Test GPU vector.

      dca::linalg::matrixop::multiplyDiagonalRight(da, dd, db);
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(db);

      for (int j = 0; j < a.nrCols(); ++j)
        for (int i = 0; i < a.nrRows(); ++i) {
          EXPECT_GE(10 * this->epsilon, std::abs(b(i, j) - d[j] * a(i, j)));
        }
    }
  }
}

TYPED_TEST(MatrixopComplexGPUTest, Trsm) {
  using ScalarType = TypeParam;
  auto val_a = [](int i, int j) { return ScalarType(3 * i - 2 * j, 1 - i * i + j); };
  auto val_b = [](int i, int j) { return ScalarType(i - 2 * j * j, 4 * i - 3 * j); };
  std::pair<int, int> size_a(3, 3);
  std::pair<int, int> size_b(3, 4);

  char uplo_opts[] = {'U', 'L'};
  char diag_opts[] = {'U', 'N'};

  for (auto uplo : uplo_opts) {
    for (auto diag : diag_opts) {
      int m = size_b.first;
      int n = size_b.second;
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size_b);

      testing::setMatrixElements(a, val_a);
      testing::setMatrixElements(b, val_b);

      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> da(a);
      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> db(b);

      dca::linalg::matrixop::trsm(uplo, diag, da, db);

      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b_test(db);

      dca::linalg::blas::trmm("L", &uplo, "N", &diag, m, n, ScalarType(1), a.ptr(),
                              a.leadingDimension(), b_test.ptr(), b_test.leadingDimension());

      for (int j = 0; j < b.nrCols(); ++j)
        for (int i = 0; i < b.nrRows(); ++i) {
          EXPECT_GE(500 * this->epsilon, std::abs(b(i, j) - b_test(i, j)));
        }
    }
  }
}

TYPED_TEST(MatrixopComplexGPUTest, Laset) {
  using ScalarType = TypeParam;
  std::pair<int, int> size(41, 35);
  ScalarType diag(3.4, 1.11);
  ScalarType offdiag(-1.4, 0.1);

  dca::linalg::Matrix<ScalarType, dca::linalg::GPU> da(size);
  dca::linalg::lapack::UseDevice<dca::linalg::GPU>::laset(size.first, size.second, offdiag, diag,
                                                          da.ptr(), da.leadingDimension(), 0, 0);

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(da);

  for (int j = 0; j < a.nrCols(); ++j)
    for (int i = 0; i < a.nrRows(); ++i) {
      if (i == j)
        EXPECT_EQ(diag, a(i, j));
      else
        EXPECT_EQ(offdiag, a(i, j));
    }
}

TYPED_TEST(MatrixopComplexGPUTest, Inverse) {
  dca::linalg::util::initializeMagma();

  using ScalarType = TypeParam;
  int size = 6;
  auto val = [size](int i, int j) {
    return std::polar<typename ScalarType::value_type>(2. / ((4 + size + i - j) % size + 1), i + j);
  };
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size);
  testing::setMatrixElements(mat, val);

  dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dinvmat(mat);
  dca::linalg::matrixop::inverse(dinvmat);
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> invmat(dinvmat);

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> res(size);
  dca::linalg::matrixop::gemm(mat, invmat, res);

  for (int j = 0; j < mat.nrCols(); ++j) {
    for (int i = 0; i < mat.nrRows(); ++i) {
      if (i == j)
        EXPECT_GE(500 * this->epsilon, std::abs(ScalarType(1) - res(i, j)));
      else
        EXPECT_GE(500 * this->epsilon, std::abs(res(i, j)));
    }
  }
}

TYPED_TEST(MatrixopComplexGPUTest, RemoveRowCol) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2(41, 35);
  auto val = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size2);
  testing::setMatrixElements(mat, val);

  for (int ii : {0, 1, size2.first - 1}) {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dmat(mat);

    dca::linalg::matrixop::removeRow(dmat, ii);
    EXPECT_EQ(mat.nrRows() - 1, dmat.nrRows());
    EXPECT_EQ(mat.nrCols(), dmat.nrCols());

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat_test(dmat);

    for (int j = 0; j < mat.nrCols(); ++j) {
      for (int i = 0; i < ii; ++i)
        EXPECT_EQ(mat(i, j), mat_test(i, j));
      for (int i = ii + 1; i < mat.nrRows(); ++i)
        EXPECT_EQ(mat(i, j), mat_test(i - 1, j));
    }
  }
  for (int jj : {0, 1, size2.second - 1}) {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dmat(mat);

    dca::linalg::matrixop::removeCol(dmat, jj);
    EXPECT_EQ(mat.nrRows(), dmat.nrRows());
    EXPECT_EQ(mat.nrCols() - 1, dmat.nrCols());

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat_test(dmat);

    for (int i = 0; i < mat.nrRows(); ++i) {
      for (int j = 0; j < jj; ++j)
        EXPECT_EQ(mat(i, j), mat_test(i, j));
      for (int j = jj + 1; j < mat.nrCols(); ++j)
        EXPECT_EQ(mat(i, j), mat_test(i, j - 1));
    }
  }
  for (int ii : {0, 1, size2.first - 1}) {
    for (int jj : {0, 1, size2.second - 1}) {
      dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dmat(mat);

      dca::linalg::matrixop::removeRowAndCol(dmat, ii, jj);
      EXPECT_EQ(mat.nrRows() - 1, dmat.nrRows());
      EXPECT_EQ(mat.nrCols() - 1, dmat.nrCols());

      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat_test(dmat);

      for (int j = 0; j < jj; ++j) {
        for (int i = 0; i < ii; ++i)
          EXPECT_EQ(mat(i, j), mat_test(i, j));
        for (int i = ii + 1; i < mat.nrRows(); ++i)
          EXPECT_EQ(mat(i, j), mat_test(i - 1, j));
      }
      for (int j = jj + 1; j < mat.nrCols(); ++j) {
        for (int i = 0; i < ii; ++i)
          EXPECT_EQ(mat(i, j), mat_test(i, j - 1));
        for (int i = ii + 1; i < mat.nrRows(); ++i)
          EXPECT_EQ(mat(i, j), mat_test(i - 1, j - 1));
      }
    }
  }
  for (int ii : {0, 1, std::min(size2.first, size2.second) - 1}) {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dmat(mat);

    dca::linalg::matrixop::removeRowAndCol(dmat, ii);
    EXPECT_EQ(mat.nrRows() - 1, dmat.nrRows());
    EXPECT_EQ(mat.nrCols() - 1, dmat.nrCols());

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat_test(dmat);

    for (int j = 0; j < ii; ++j) {
      for (int i = 0; i < ii; ++i)
        EXPECT_EQ(mat(i, j), mat_test(i, j));
      for (int i = ii + 1; i < mat.nrRows(); ++i)
        EXPECT_EQ(mat(i, j), mat_test(i - 1, j));
    }
    for (int j = ii + 1; j < mat.nrCols(); ++j) {
      for (int i = 0; i < ii; ++i)
        EXPECT_EQ(mat(i, j), mat_test(i, j - 1));
      for (int i = ii + 1; i < mat.nrRows(); ++i)
        EXPECT_EQ(mat(i, j), mat_test(i - 1, j - 1));
    }
  }
}

TYPED_TEST(MatrixopComplexGPUTest, CopyRow) {
  using ScalarType = TypeParam;
  const ScalarType checked(-1000);
  std::pair<int, int> size2_a(41, 35);
  std::pair<int, int> size2_b(46, 35);

  dca::linalg::Vector<int, dca::linalg::CPU> i_sources(3);
  i_sources[0] = 0;
  i_sources[1] = 2;
  i_sources[2] = 3;
  dca::linalg::Vector<int, dca::linalg::GPU> di_sources(i_sources);

  dca::linalg::Vector<int, dca::linalg::CPU> i_dests(4);
  i_dests[0] = 2;
  i_dests[1] = 5;
  i_dests[2] = 3;
  i_dests[3] = 4;
  dca::linalg::Vector<int, dca::linalg::GPU> di_dests(i_dests);

  auto val_a = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };
  auto val_b = [](int i, int j) { return 10 * i + j + 100; };
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size2_b);
  testing::setMatrixElements(b, val_b);

  dca::linalg::Matrix<ScalarType, dca::linalg::GPU> da(a);
  {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(b);

    for (int i = 0; i < i_sources.size(); ++i)
      dca::linalg::matrixop::copyRow(da, i_sources[i], dc, i_dests[i]);

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(dc);

    for (int j = 0; j < a.nrCols(); ++j) {
      for (int i = 0; i < i_sources.size(); ++i) {
        EXPECT_EQ(a(i_sources[i], j), c(i_dests[i], j));
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(i_dests[i], j) = checked;
      }
      for (int i = 0; i < b.nrRows(); ++i)
        if (c(i, j) != checked)
          EXPECT_EQ(b(i, j), c(i, j));
    }
  }
  {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(b);

    dca::linalg::matrixop::copyRows(da, di_sources, dc, di_dests);

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(dc);

    for (int j = 0; j < a.nrCols(); ++j) {
      for (int i = 0; i < i_sources.size(); ++i) {
        EXPECT_EQ(a(i_sources[i], j), c(i_dests[i], j));
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(i_dests[i], j) = checked;
      }
      for (int i = 0; i < b.nrRows(); ++i)
        if (c(i, j) != checked)
          EXPECT_EQ(b(i, j), c(i, j));
    }
  }
}

TYPED_TEST(MatrixopComplexGPUTest, CopyCol) {
  using ScalarType = TypeParam;
  const ScalarType checked(-1000);
  std::pair<int, int> size2_a(35, 34);
  std::pair<int, int> size2_b(35, 36);

  dca::linalg::Vector<int, dca::linalg::CPU> j_sources(3);
  j_sources[0] = 0;
  j_sources[1] = 2;
  j_sources[2] = 3;
  dca::linalg::Vector<int, dca::linalg::GPU> dj_sources(j_sources);

  dca::linalg::Vector<int, dca::linalg::CPU> j_dests(4);
  j_dests[0] = 2;
  j_dests[1] = 5;
  j_dests[2] = 3;
  j_dests[3] = 4;
  dca::linalg::Vector<int, dca::linalg::GPU> dj_dests(j_dests);

  auto val_a = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };
  auto val_b = [](int i, int j) { return 10 * i + j + 100; };
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size2_b);
  testing::setMatrixElements(b, val_b);

  dca::linalg::Matrix<ScalarType, dca::linalg::GPU> da(a);

  {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(b);

    for (int j = 0; j < j_sources.size(); ++j)
      dca::linalg::matrixop::copyCol(da, j_sources[j], dc, j_dests[j]);

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(dc);

    for (int i = 0; i < a.nrRows(); ++i) {
      for (int j = 0; j < j_sources.size(); ++j) {
        EXPECT_EQ(a(i, j_sources[j]), c(i, j_dests[j]));
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(i, j_dests[j]) = checked;
      }
      for (int j = 0; j < b.nrCols(); ++j)
        if (c(i, j) != checked)
          EXPECT_EQ(b(i, j), c(i, j));
    }
  }
  {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(b);

    dca::linalg::matrixop::copyCols(da, dj_sources, dc, dj_dests);

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(dc);

    for (int i = 0; i < a.nrRows(); ++i) {
      for (int j = 0; j < j_sources.size(); ++j) {
        EXPECT_EQ(a(i, j_sources[j]), c(i, j_dests[j]));
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(i, j_dests[j]) = checked;
      }
      for (int j = 0; j < b.nrCols(); ++j)
        if (c(i, j) != checked)
          EXPECT_EQ(b(i, j), c(i, j));
    }
  }
}

TYPED_TEST(MatrixopComplexGPUTest, ScaleRow) {
  using ScalarType = TypeParam;
  const ScalarType checked(-1000);
  std::pair<int, int> size2_a(41, 35);

  dca::linalg::Vector<int, dca::linalg::CPU> is(3);
  is[0] = 0;
  is[1] = 2;
  is[2] = 3;
  dca::linalg::Vector<int, dca::linalg::GPU> dis(is);

  dca::linalg::Vector<ScalarType, dca::linalg::CPU> vals(3);
  vals[0] = 3.4;
  vals[1] = -1.2;
  vals[2] = 7.7;
  dca::linalg::Vector<ScalarType, dca::linalg::GPU> dvals(vals);

  auto val_a = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);

  {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(a);

    for (int i = 0; i < is.size(); ++i)
      dca::linalg::matrixop::scaleRow(dc, is[i], vals[i]);

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(dc);

    for (int j = 0; j < a.nrCols(); ++j) {
      for (int i = 0; i < is.size(); ++i) {
        EXPECT_GE(10 * this->epsilon, std::abs(vals[i] * a(is[i], j) - c(is[i], j)));
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(is[i], j) = checked;
      }
      for (int i = 0; i < a.nrRows(); ++i)
        if (c(i, j) != checked)
          EXPECT_EQ(a(i, j), c(i, j));
    }
  }
  {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(a);

    dca::linalg::matrixop::scaleRows(dc, dis, dvals);

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(dc);

    for (int j = 0; j < a.nrCols(); ++j) {
      for (int i = 0; i < is.size(); ++i) {
        EXPECT_GE(10 * this->epsilon, std::abs(vals[i] * a(is[i], j) - c(is[i], j)));
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(is[i], j) = checked;
      }
      for (int i = 0; i < a.nrRows(); ++i)
        if (c(i, j) != checked)
          EXPECT_EQ(a(i, j), c(i, j));
    }
  }
}

TYPED_TEST(MatrixopComplexGPUTest, ScaleCol) {
  using ScalarType = TypeParam;
  const ScalarType checked(-1000);
  std::pair<int, int> size2_a(31, 45);

  dca::linalg::Vector<int, dca::linalg::CPU> js(3);
  js[0] = 0;
  js[1] = 2;
  js[2] = 3;

  dca::linalg::Vector<ScalarType, dca::linalg::CPU> vals(3);
  vals[0] = 3.4;
  vals[1] = -1.2;
  vals[2] = 7.7;

  auto val_a = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);

  {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(a);

    for (int j = 0; j < js.size(); ++j)
      dca::linalg::matrixop::scaleCol(dc, js[j], vals[j]);

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(dc);

    for (int i = 0; i < a.nrRows(); ++i) {
      for (int j = 0; j < js.size(); ++j) {
        EXPECT_GE(10 * this->epsilon, std::abs(vals[j] * a(i, js[j]) - c(i, js[j])));
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(i, js[j]) = checked;
      }
      for (int j = 0; j < a.nrCols(); ++j)
        if (c(i, j) != checked)
          EXPECT_EQ(a(i, j), c(i, j));
    }
  }
}

TYPED_TEST(MatrixopComplexGPUTest, SwapRow) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2_a(41, 35);

  dca::linalg::Vector<int, dca::linalg::CPU> i_1(3);
  i_1[0] = 0;
  i_1[1] = 2;
  i_1[2] = 3;
  dca::linalg::Vector<int, dca::linalg::GPU> di_1(i_1);

  dca::linalg::Vector<int, dca::linalg::CPU> i_2(3);
  i_2[0] = 1;
  i_2[1] = 5;
  i_2[2] = 4;
  dca::linalg::Vector<int, dca::linalg::GPU> di_2(i_2);

  auto val_a = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);

  {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(a);

    for (int i = 0; i < i_1.size(); ++i) {
      dca::linalg::matrixop::swapRow(a, i_1[i], i_2[i]);
      dca::linalg::matrixop::swapRow(dc, i_1[i], i_2[i]);
    }

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(dc);

    for (int j = 0; j < a.nrCols(); ++j)
      for (int i = 0; i < a.nrRows(); ++i)
        EXPECT_EQ(a(i, j), c(i, j));
  }
  {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(a);

    for (int i = 0; i < i_1.size(); ++i)
      dca::linalg::matrixop::swapRow(a, i_1[i], i_2[i]);

    dca::linalg::matrixop::swapRows(dc, di_1, di_2);

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(dc);

    for (int j = 0; j < a.nrCols(); ++j)
      for (int i = 0; i < a.nrRows(); ++i)
        EXPECT_EQ(a(i, j), c(i, j));
  }
}

TYPED_TEST(MatrixopComplexGPUTest, SwapCol) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2_a(41, 35);

  dca::linalg::Vector<int, dca::linalg::CPU> j_1(3);
  j_1[0] = 0;
  j_1[1] = 2;
  j_1[2] = 3;
  dca::linalg::Vector<int, dca::linalg::GPU> dj_1(j_1);

  dca::linalg::Vector<int, dca::linalg::CPU> j_2(3);
  j_2[0] = 1;
  j_2[1] = 5;
  j_2[2] = 4;
  dca::linalg::Vector<int, dca::linalg::GPU> dj_2(j_2);

  auto val_a = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);

  {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(a);

    for (int j = 0; j < j_1.size(); ++j) {
      dca::linalg::matrixop::swapCol(a, j_1[j], j_2[j]);
      dca::linalg::matrixop::swapCol(dc, j_1[j], j_2[j]);
    }

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(dc);

    for (int j = 0; j < a.nrCols(); ++j)
      for (int i = 0; i < a.nrRows(); ++i)
        EXPECT_EQ(a(i, j), c(i, j));
  }
  {
    dca::linalg::Matrix<ScalarType, dca::linalg::GPU> dc(a);

    for (int j = 0; j < j_1.size(); ++j)
      dca::linalg::matrixop::swapCol(a, j_1[j], j_2[j]);

    dca::linalg::matrixop::swapCols(dc, dj_1, dj_2);

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(dc);

    for (int j = 0; j < a.nrCols(); ++j)
      for (int i = 0; i < a.nrRows(); ++i)
        EXPECT_EQ(a(i, j), c(i, j));
  }
}
