// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the matrix operations with the Matrix<CPU> class with real element type.

#include "dca/linalg/matrixop.hpp"
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include "gtest/gtest.h"
#include "dca/linalg/lapack/use_device.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/blas/blas3.hpp"
#include "dca/math/random/random.hpp"
#include "cpu_test_util.hpp"
#include "matrixop_reference.hpp"

template <typename ScalarType>
class MatrixopRealCPUTest : public ::testing::Test {
public:
  static const ScalarType epsilon;
};
template <typename ScalarType>
const ScalarType MatrixopRealCPUTest<ScalarType>::epsilon = std::numeric_limits<ScalarType>::epsilon();

typedef ::testing::Types<float, double> FloatingPointTypes;
TYPED_TEST_CASE(MatrixopRealCPUTest, FloatingPointTypes);

TYPED_TEST(MatrixopRealCPUTest, Gemv) {
  using ScalarType = TypeParam;
  auto val_a = [](int i, int j) { return 3 * i - 2 * j; };
  auto val_x = [](int i) { return 4 * i; };
  auto val_y = [](int i) { return 2 * i * i; };

  int m = 2;
  int n = 3;

  char trans_opts[] = {'N', 'T', 'C'};

  for (auto transa : trans_opts) {
    std::pair<int, int> size_a;
    if (transa == 'N') {
      size_a.first = m;
      size_a.second = n;
    }
    else {
      size_a.first = n;
      size_a.second = m;
    }

    {
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
      dca::linalg::Vector<ScalarType, dca::linalg::CPU> x(n);
      dca::linalg::Vector<ScalarType, dca::linalg::CPU> y(m);

      testing::setMatrixElements(a, val_a);
      testing::setVectorElements(x, val_x);
      testing::setVectorElements(y, val_y);

      dca::linalg::Vector<ScalarType, dca::linalg::CPU> yref(y);

      testing::refGemv(transa, ScalarType(1), a, x, ScalarType(0), yref);
      dca::linalg::matrixop::gemv(transa, a, x, y);

      for (int i = 0; i < y.size(); ++i) {
        EXPECT_NEAR(yref[i], y[i], 500 * this->epsilon);
      }
    }
    {
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
      dca::linalg::Vector<ScalarType, dca::linalg::CPU> x(n);
      dca::linalg::Vector<ScalarType, dca::linalg::CPU> y(m);

      testing::setMatrixElements(a, val_a);
      testing::setVectorElements(x, val_x);
      testing::setVectorElements(y, val_y);

      dca::linalg::Vector<ScalarType, dca::linalg::CPU> yref(y);

      ScalarType alpha = .57;
      ScalarType beta = -.712;

      testing::refGemv(transa, alpha, a, x, beta, yref);
      dca::linalg::matrixop::gemv(transa, alpha, a, x, beta, y);

      for (int i = 0; i < y.size(); ++i) {
        EXPECT_NEAR(yref[i], y[i], 500 * this->epsilon);
      }
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, Gemm) {
  using ScalarType = TypeParam;
  auto val_a = [](int i, int j) { return 3 * i - 2 * j; };
  auto val_b = [](int i, int j) { return 4 * i - 3 * j; };
  auto val_c = [](int i, int j) { return 2 * i - j; };
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

      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> cref(c);

      testing::refGemm('N', 'N', ScalarType(1), a, b, ScalarType(0), cref);
      dca::linalg::matrixop::gemm(a, b, c);

      for (int j = 0; j < c.nrCols(); ++j)
        for (int i = 0; i < c.nrRows(); ++i) {
          EXPECT_NEAR(cref(i, j), c(i, j), 500 * this->epsilon);
        }
    }
    {
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size_b);
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(size_c);

      testing::setMatrixElements(a, val_a);
      testing::setMatrixElements(b, val_b);
      testing::setMatrixElements(c, val_c);

      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> cref(c);

      ScalarType alpha = .57;
      ScalarType beta = -.712;

      testing::refGemm('N', 'N', alpha, a, b, beta, cref);
      dca::linalg::matrixop::gemm(alpha, a, b, beta, c);

      for (int j = 0; j < c.nrCols(); ++j)
        for (int i = 0; i < c.nrRows(); ++i) {
          EXPECT_NEAR(cref(i, j), c(i, j), 500 * this->epsilon);
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

          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> cref(c);

          testing::refGemm(transa, transb, ScalarType(1), a, b, ScalarType(0), cref);
          dca::linalg::matrixop::gemm(transa, transb, a, b, c);

          for (int j = 0; j < c.nrCols(); ++j)
            for (int i = 0; i < c.nrRows(); ++i) {
              EXPECT_NEAR(cref(i, j), c(i, j), 500 * this->epsilon);
            }
        }
        {
          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size_b);
          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(size_c);

          testing::setMatrixElements(a, val_a);
          testing::setMatrixElements(b, val_b);
          testing::setMatrixElements(c, val_c);

          dca::linalg::Matrix<ScalarType, dca::linalg::CPU> cref(c);

          ScalarType alpha = .57;
          ScalarType beta = -.712;

          testing::refGemm(transa, transb, alpha, a, b, beta, cref);
          dca::linalg::matrixop::gemm(transa, transb, alpha, a, b, beta, c);

          for (int j = 0; j < c.nrCols(); ++j)
            for (int i = 0; i < c.nrRows(); ++i) {
              EXPECT_NEAR(cref(i, j), c(i, j), 500 * this->epsilon);
            }
        }
      }
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, MultiplyDiagonal) {
  using ScalarType = TypeParam;
  std::pair<int, int> size_a(3, 5);
  auto val_a = [](int i, int j) { return 3 * i - 2 * j; };
  auto val_d = [](int i) { return 1 - i; };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
  testing::setMatrixElements(a, val_a);
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size_a);
  {
    dca::linalg::Vector<ScalarType, dca::linalg::CPU> d(a.nrRows());
    testing::setVectorElements(d, val_d);

    dca::linalg::matrixop::multiplyDiagonalLeft(d, a, b);

    for (int j = 0; j < a.nrCols(); ++j)
      for (int i = 0; i < a.nrRows(); ++i) {
        EXPECT_NEAR(d[i] * a(i, j), b(i, j), 10 * this->epsilon);
      }
  }
  {
    dca::linalg::Vector<ScalarType, dca::linalg::CPU> d(a.nrCols());
    testing::setVectorElements(d, val_d);

    dca::linalg::matrixop::multiplyDiagonalRight(a, d, b);

    for (int j = 0; j < a.nrCols(); ++j)
      for (int i = 0; i < a.nrRows(); ++i) {
        EXPECT_NEAR(d[j] * a(i, j), b(i, j), 10 * this->epsilon);
      }
  }
}

TYPED_TEST(MatrixopRealCPUTest, Trsm) {
  using ScalarType = TypeParam;
  auto val_a = [](int i, int j) { return 1 + 3 * i - 2 * j; };
  auto val_b = [](int i, int j) { return 4 * i - 3 * j; };
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

      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b_orig(b);

      dca::linalg::matrixop::trsm(uplo, diag, a, b);
      dca::linalg::blas::trmm("L", &uplo, "N", &diag, m, n, ScalarType(1), a.ptr(),
                              a.leadingDimension(), b.ptr(), b.leadingDimension());

      for (int j = 0; j < b.nrCols(); ++j)
        for (int i = 0; i < b.nrRows(); ++i) {
          EXPECT_NEAR(b_orig(i, j), b(i, j), 500 * this->epsilon);
        }
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, Eigensolver) {
  using ScalarType = TypeParam;
  int size = 2;
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size);
  mat(0, 0) = 2.;
  mat(0, 1) = 0.;
  mat(1, 0) = 1.;
  mat(1, 1) = 1.;

  dca::linalg::Vector<ScalarType, dca::linalg::CPU> wr;
  dca::linalg::Vector<ScalarType, dca::linalg::CPU> wi;
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> vl;
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> vr;

  dca::linalg::matrixop::eigensolver('V', 'V', mat, wr, wi, vl, vr);
  EXPECT_EQ(wr.size(), size);
  EXPECT_EQ(wi.size(), size);
  EXPECT_EQ(vl.size(), std::make_pair(size, size));
  EXPECT_EQ(vr.size(), std::make_pair(size, size));

  EXPECT_NEAR(1., wr[0], 100 * this->epsilon);
  EXPECT_NEAR(0., wi[0], 100 * this->epsilon);
  EXPECT_NEAR(-1, vl(0, 0) / vl(1, 0), 100 * this->epsilon);
  EXPECT_NEAR(0., vr(0, 0), 100 * this->epsilon);
  EXPECT_NEAR(1., vr(1, 0), 100 * this->epsilon);

  EXPECT_NEAR(2., wr[1], 100 * this->epsilon);
  EXPECT_NEAR(0., wi[1], 100 * this->epsilon);
  EXPECT_NEAR(1., vl(0, 1), 100 * this->epsilon);
  EXPECT_NEAR(0., vl(1, 1), 100 * this->epsilon);
  EXPECT_NEAR(1., vr(0, 1) / vr(1, 1), 100 * this->epsilon);
}

TYPED_TEST(MatrixopRealCPUTest, EigensolverSymmetric) {
  using ScalarType = TypeParam;
  int size = 2;
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size);
  mat(0, 0) = 2.;
  mat(0, 1) = 1.;
  mat(1, 0) = 1.;
  mat(1, 1) = 2.;

  dca::linalg::Vector<ScalarType, dca::linalg::CPU> w;
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> v;

  dca::linalg::matrixop::eigensolverSymmetric('V', 'U', mat, w, v);
  EXPECT_EQ(w.size(), size);
  EXPECT_EQ(v.size(), std::make_pair(size, size));

  EXPECT_NEAR(1., w[0], 100 * this->epsilon);
  EXPECT_NEAR(-1., v(0, 0) / v(1, 0), 100 * this->epsilon);

  EXPECT_NEAR(3., w[1], 100 * this->epsilon);
  EXPECT_NEAR(1., v(0, 1) / v(1, 1), 100 * this->epsilon);

  dca::linalg::matrixop::eigensolverHermitian('V', 'U', mat, w, v);
  EXPECT_EQ(w.size(), size);
  EXPECT_EQ(v.size(), std::make_pair(size, size));

  EXPECT_NEAR(1., w[0], 100 * this->epsilon);
  EXPECT_NEAR(-1., v(0, 0) / v(1, 0), 100 * this->epsilon);

  EXPECT_NEAR(3., w[1], 100 * this->epsilon);
  EXPECT_NEAR(1., v(0, 1) / v(1, 1), 100 * this->epsilon);
}

TYPED_TEST(MatrixopRealCPUTest, Laset) {
  using ScalarType = TypeParam;
  std::pair<int, int> size(3, 5);
  ScalarType diag = 3.4;
  ScalarType offdiag = -1.4;

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size);
  dca::linalg::lapack::UseDevice<dca::linalg::CPU>::laset(size.first, size.second, offdiag, diag,
                                                          a.ptr(), a.leadingDimension(), -1, -1);

  for (int j = 0; j < a.nrCols(); ++j)
    for (int i = 0; i < a.nrRows(); ++i) {
      if (i == j)
        EXPECT_EQ(diag, a(i, j));
      else
        EXPECT_EQ(offdiag, a(i, j));
    }
}

TYPED_TEST(MatrixopRealCPUTest, InsertRowCol) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2(31, 31);
  auto val = [](int i, int j) { return 10 * i + j; };
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size2);
  testing::setMatrixElements(mat, val);

  for (int ii : {0, 1, size2.first}) {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat_test(mat);

    dca::linalg::matrixop::insertRow(mat_test, ii);
    for (int j = 0; j < mat.nrCols(); ++j) {
      for (int i = 0; i < ii; ++i)
        EXPECT_EQ(mat(i, j), mat_test(i, j));
      for (int i = ii; i < mat.nrRows(); ++i)
        EXPECT_EQ(mat(i, j), mat_test(i + 1, j));

      EXPECT_EQ(0, mat_test(ii, j));
    }
  }
  for (int jj : {0, 1, size2.second}) {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat_test(mat);

    dca::linalg::matrixop::insertCol(mat_test, jj);
    for (int i = 0; i < mat.nrRows(); ++i) {
      for (int j = 0; j < jj; ++j)
        EXPECT_EQ(mat(i, j), mat_test(i, j));
      for (int j = jj; j < mat.nrCols(); ++j)
        EXPECT_EQ(mat(i, j), mat_test(i, j + 1));

      EXPECT_EQ(0, mat_test(i, jj));
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, Inverse) {
  using ScalarType = TypeParam;
  int size = 6;
  auto val = [size](int i, int j) { return ScalarType(2. / ((4 + size + i - j) % size + 1)); };
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size);
  testing::setMatrixElements(mat, val);

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> invmat(mat);
  dca::linalg::matrixop::inverse(invmat);

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> res(size);
  dca::linalg::matrixop::gemm(mat, invmat, res);

  for (int j = 0; j < mat.nrCols(); ++j) {
    for (int i = 0; i < mat.nrRows(); ++i) {
      if (i == j)
        EXPECT_NEAR(1, res(i, j), 500 * this->epsilon);
      else
        EXPECT_NEAR(0, res(i, j), 500 * this->epsilon);
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, PseudoInverse) {
  using ScalarType = TypeParam;
  auto val = [](int i, int j) { return ScalarType(2. / ((5 + i - j) % 3 + 1)); };
  auto val0 = [](int, int) { return 0; };
  {
    std::pair<int, int> size(2, 3);
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size);
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> invmat;
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> res(std::min(size.first, size.second));

    testing::setMatrixElements(mat, val);
    dca::linalg::matrixop::pseudoInverse(mat, invmat);
    if (size.first <= size.second)
      dca::linalg::matrixop::gemm(mat, invmat, res);
    else
      dca::linalg::matrixop::gemm(invmat, mat, res);

    for (int j = 0; j < res.nrCols(); ++j) {
      for (int i = 0; i < res.nrRows(); ++i) {
        if (i == j)
          EXPECT_NEAR(1, res(i, j), 500 * this->epsilon);
        else
          EXPECT_NEAR(0, res(i, j), 500 * this->epsilon);
      }
    }

    // Check eigenvalue exclusion below threshold.
    testing::setMatrixElements(mat, val0);
    mat(0, 0) = 20.;
    mat(1, 1) = .1;
    dca::linalg::matrixop::pseudoInverse(mat, invmat, 1e-6);
    if (size.first <= size.second)
      dca::linalg::matrixop::gemm(mat, invmat, res);
    else
      dca::linalg::matrixop::gemm(invmat, mat, res);

    for (int j = 0; j < res.nrCols(); ++j) {
      for (int i = 0; i < res.nrRows(); ++i) {
        if (i == j)
          EXPECT_NEAR(1, res(i, j), 100 * this->epsilon);
        else
          EXPECT_NEAR(0, res(i, j), 100 * this->epsilon);
      }
    }

    // Check eigenvalue exclusion above threshold.
    dca::linalg::matrixop::pseudoInverse(mat, invmat, 3.9e-4);
    if (size.first <= size.second)
      dca::linalg::matrixop::gemm(mat, invmat, res);
    else
      dca::linalg::matrixop::gemm(invmat, mat, res);

    for (int j = 0; j < res.nrCols(); ++j) {
      for (int i = 0; i < res.nrRows(); ++i) {
        if (i == j && i == 0)
          EXPECT_NEAR(1, res(i, j), 100 * this->epsilon);
        else
          EXPECT_NEAR(0, res(i, j), 100 * this->epsilon);
      }
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, RemoveRowCol) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2(3, 4);
  auto val = [](int i, int j) { return 10 * i + j; };
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size2);
  testing::setMatrixElements(mat, val);

  for (int ii : {0, 1, size2.first - 1}) {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat_test(mat);

    dca::linalg::matrixop::removeRow(mat_test, ii);
    EXPECT_EQ(mat.nrRows() - 1, mat_test.nrRows());
    EXPECT_EQ(mat.nrCols(), mat_test.nrCols());

    for (int j = 0; j < mat.nrCols(); ++j) {
      for (int i = 0; i < ii; ++i)
        EXPECT_EQ(mat(i, j), mat_test(i, j));
      for (int i = ii + 1; i < mat.nrRows(); ++i)
        EXPECT_EQ(mat(i, j), mat_test(i - 1, j));
    }
  }
  for (int jj : {0, 1, size2.second - 1}) {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat_test(mat);

    dca::linalg::matrixop::removeCol(mat_test, jj);
    EXPECT_EQ(mat.nrRows(), mat_test.nrRows());
    EXPECT_EQ(mat.nrCols() - 1, mat_test.nrCols());

    for (int i = 0; i < mat.nrRows(); ++i) {
      for (int j = 0; j < jj; ++j)
        EXPECT_EQ(mat(i, j), mat_test(i, j));
      for (int j = jj + 1; j < mat.nrCols(); ++j)
        EXPECT_EQ(mat(i, j), mat_test(i, j - 1));
    }
  }
  for (int ii : {0, 1, size2.first - 1}) {
    for (int jj : {0, 1, size2.second - 1}) {
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat_test(mat);

      dca::linalg::matrixop::removeRowAndCol(mat_test, ii, jj);
      EXPECT_EQ(mat.nrRows() - 1, mat_test.nrRows());
      EXPECT_EQ(mat.nrCols() - 1, mat_test.nrCols());

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
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat_test(mat);

    dca::linalg::matrixop::removeRowAndCol(mat_test, ii);
    EXPECT_EQ(mat.nrRows() - 1, mat_test.nrRows());
    EXPECT_EQ(mat.nrCols() - 1, mat_test.nrCols());

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

TEST(MatrixopRealCPUTest, RemoveCols) {
  dca::linalg::Matrix<int, dca::linalg::CPU> mat(4);
  dca::linalg::Matrix<int, dca::linalg::CPU> mat_test(std::make_pair(4, 2));

  auto fill_func = [](int i, int j) { return 10 * i + j; };
  testing::setMatrixElements(mat, fill_func);

  mat_test(0, 0) = 0, mat_test(0, 1) = 3;
  mat_test(1, 0) = 10, mat_test(1, 1) = 13;
  mat_test(2, 0) = 20, mat_test(2, 1) = 23;
  mat_test(3, 0) = 30, mat_test(3, 1) = 33;

  dca::linalg::matrixop::removeCols(mat, 1, 2);
  EXPECT_EQ(mat, mat_test);

  dca::linalg::Matrix<int, dca::linalg::CPU> null_matrix(0);
  dca::linalg::matrixop::removeCols(mat, 0, mat.nrCols() - 1);
  EXPECT_EQ(null_matrix, mat);
}

TEST(MatrixopRealCPUTest, RemoveRows) {
  dca::linalg::Matrix<int, dca::linalg::CPU> mat(4);
  dca::linalg::Matrix<int, dca::linalg::CPU> mat_test(std::make_pair(2, 4));

  auto fill_func = [](int i, int j) { return 10 * i + j; };
  testing::setMatrixElements(mat, fill_func);

  mat_test(0, 0) = 0, mat_test(0, 1) = 1, mat_test(0, 2) = 2, mat_test(0, 3) = 3;
  mat_test(1, 0) = 30, mat_test(1, 1) = 31, mat_test(1, 2) = 32, mat_test(1, 3) = 33;

  dca::linalg::matrixop::removeRows(mat, 1, 2);

  EXPECT_EQ(mat, mat_test);

  dca::linalg::Matrix<int, dca::linalg::CPU> null_matrix(0);
  dca::linalg::matrixop::removeRows(mat, 0, mat.nrRows() - 1);
  EXPECT_EQ(mat, null_matrix);
}

TEST(MatrixopRealCPUTest, RemoveRowsandColumns) {
  dca::linalg::Matrix<int, dca::linalg::CPU> mat(4);
  dca::linalg::Matrix<int, dca::linalg::CPU> mat_test(2);

  auto fill_func = [](int i, int j) { return 10 * i + j; };
  testing::setMatrixElements(mat, fill_func);

  mat_test(0, 0) = 0, mat_test(0, 1) = 1;
  mat_test(1, 0) = 10, mat_test(1, 1) = 11;

  dca::linalg::matrixop::removeRowsAndCols(mat, 2, 3);

  EXPECT_EQ(mat, mat_test);
}

TYPED_TEST(MatrixopRealCPUTest, CopyRow) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2_a(4, 3);
  std::pair<int, int> size2_b(6, 3);

  dca::linalg::Vector<int, dca::linalg::CPU> i_sources(3);
  i_sources[0] = 0;
  i_sources[1] = 2;
  i_sources[2] = 3;

  dca::linalg::Vector<int, dca::linalg::CPU> i_dests(4);
  i_dests[0] = 2;
  i_dests[1] = 5;
  i_dests[2] = 3;
  i_dests[3] = 4;

  auto val_a = [](int i, int j) { return 10 * i + j; };
  auto val_b = [](int i, int j) { return 10 * i + j + 100; };
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size2_b);
  testing::setMatrixElements(b, val_b);

  {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(b);

    for (int i = 0; i < i_sources.size(); ++i)
      dca::linalg::matrixop::copyRow(a, i_sources[i], c, i_dests[i]);

    for (int j = 0; j < a.nrCols(); ++j) {
      for (int i = 0; i < i_sources.size(); ++i) {
        EXPECT_EQ(a(i_sources[i], j), c(i_dests[i], j));
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(i_dests[i], j) = -1000;
      }
      for (int i = 0; i < b.nrRows(); ++i)
        if (c(i, j) != -1000)
          EXPECT_EQ(b(i, j), c(i, j));
    }
  }
  {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(b);

    dca::linalg::matrixop::copyRows(a, i_sources, c, i_dests);

    for (int j = 0; j < a.nrCols(); ++j) {
      for (int i = 0; i < i_sources.size(); ++i) {
        EXPECT_EQ(a(i_sources[i], j), c(i_dests[i], j));
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(i_dests[i], j) = -1000;
      }
      for (int i = 0; i < b.nrRows(); ++i)
        if (c(i, j) != -1000)
          EXPECT_EQ(b(i, j), c(i, j));
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, CopyCol) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2_a(3, 4);
  std::pair<int, int> size2_b(3, 6);

  dca::linalg::Vector<int, dca::linalg::CPU> j_sources(3);
  j_sources[0] = 0;
  j_sources[1] = 2;
  j_sources[2] = 3;

  dca::linalg::Vector<int, dca::linalg::CPU> j_dests(4);
  j_dests[0] = 2;
  j_dests[1] = 5;
  j_dests[2] = 3;
  j_dests[3] = 4;

  auto val_a = [](int i, int j) { return 10 * i + j; };
  auto val_b = [](int i, int j) { return 10 * i + j + 100; };
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size2_b);
  testing::setMatrixElements(b, val_b);

  {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(b);

    for (int j = 0; j < j_sources.size(); ++j)
      dca::linalg::matrixop::copyCol(a, j_sources[j], c, j_dests[j]);

    for (int i = 0; i < a.nrRows(); ++i) {
      for (int j = 0; j < j_sources.size(); ++j) {
        EXPECT_EQ(a(i, j_sources[j]), c(i, j_dests[j]));
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(i, j_dests[j]) = -1000;
      }
      for (int j = 0; j < b.nrCols(); ++j)
        if (c(i, j) != -1000)
          EXPECT_EQ(b(i, j), c(i, j));
    }
  }
  {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(b);

    dca::linalg::matrixop::copyCols(a, j_sources, c, j_dests);

    for (int i = 0; i < a.nrRows(); ++i) {
      for (int j = 0; j < j_sources.size(); ++j) {
        EXPECT_EQ(a(i, j_sources[j]), c(i, j_dests[j]));
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(i, j_dests[j]) = -1000;
      }
      for (int j = 0; j < b.nrCols(); ++j)
        if (c(i, j) != -1000)
          EXPECT_EQ(b(i, j), c(i, j));
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, ScaleRow) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2_a(4, 3);

  dca::linalg::Vector<int, dca::linalg::CPU> is(3);
  is[0] = 0;
  is[1] = 2;
  is[2] = 3;

  dca::linalg::Vector<ScalarType, dca::linalg::CPU> vals(3);
  vals[0] = 3.4;
  vals[1] = -1.2;
  vals[2] = 7.7;

  auto val_a = [](int i, int j) { return 10 * i + j; };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);

  {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(a);

    for (int i = 0; i < is.size(); ++i)
      dca::linalg::matrixop::scaleRow(c, is[i], vals[i]);

    for (int j = 0; j < a.nrCols(); ++j) {
      for (int i = 0; i < is.size(); ++i) {
        EXPECT_NEAR(vals[i] * a(is[i], j), c(is[i], j), 10 * this->epsilon);
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(is[i], j) = -1000;
      }
      for (int i = 0; i < a.nrRows(); ++i)
        if (c(i, j) != -1000)
          EXPECT_EQ(a(i, j), c(i, j));
    }
  }
  {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(a);

    dca::linalg::matrixop::scaleRows(c, is, vals);

    for (int j = 0; j < a.nrCols(); ++j) {
      for (int i = 0; i < is.size(); ++i) {
        EXPECT_NEAR(vals[i] * a(is[i], j), c(is[i], j), 10 * this->epsilon);
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(is[i], j) = -1000;
      }
      for (int i = 0; i < a.nrRows(); ++i)
        if (c(i, j) != -1000)
          EXPECT_EQ(a(i, j), c(i, j));
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, ScaleCol) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2_a(3, 4);

  dca::linalg::Vector<int, dca::linalg::CPU> js(3);
  js[0] = 0;
  js[1] = 2;
  js[2] = 3;

  dca::linalg::Vector<ScalarType, dca::linalg::CPU> vals(3);
  vals[0] = 3.4;
  vals[1] = -1.2;
  vals[2] = 7.7;

  auto val_a = [](int i, int j) { return 10 * i + j; };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);

  {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(a);

    for (int j = 0; j < js.size(); ++j)
      dca::linalg::matrixop::scaleCol(c, js[j], vals[j]);

    for (int i = 0; i < a.nrRows(); ++i) {
      for (int j = 0; j < js.size(); ++j) {
        EXPECT_NEAR(vals[j] * a(i, js[j]), c(i, js[j]), 10 * this->epsilon);
        // set the checked elements to -1000 to simplify the check of the unchanged elements
        c(i, js[j]) = -1000;
      }
      for (int j = 0; j < a.nrCols(); ++j)
        if (c(i, j) != -1000)
          EXPECT_EQ(a(i, j), c(i, j));
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, SwapRow) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2_a(4, 3);

  auto val_a = [](int i, int j) { return 10 * i + j; };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);

  for (int i_1 : {0, 2, 3}) {
    for (int i_2 : {1, 3}) {
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(a);

      dca::linalg::matrixop::swapRow(c, i_1, i_2);

      for (int j = 0; j < a.nrCols(); ++j) {
        for (int i = 0; i < a.nrCols(); ++i) {
          if (i == i_1) {
            EXPECT_EQ(a(i_2, j), c(i, j));
          }
          else if (i == i_2) {
            EXPECT_EQ(a(i_1, j), c(i, j));
          }
          else {
            EXPECT_EQ(a(i, j), c(i, j));
          }
        }
      }
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, SwapCol) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2_a(3, 4);

  auto val_a = [](int i, int j) { return 10 * i + j; };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);

  for (int j_1 : {0, 2, 3}) {
    for (int j_2 : {1, 3}) {
      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(a);

      dca::linalg::matrixop::swapCol(c, j_1, j_2);

      for (int i = 0; i < a.nrRows(); ++i) {
        for (int j = 0; j < a.nrCols(); ++j) {
          if (j == j_1) {
            EXPECT_EQ(a(i, j_2), c(i, j));
          }
          else if (j == j_2) {
            EXPECT_EQ(a(i, j_1), c(i, j));
          }
          else {
            EXPECT_EQ(a(i, j), c(i, j));
          }
        }
      }
    }
  }
}

TEST(MatrixopCPUTest, Difference) {
  std::pair<int, int> size2_a(5, 4);
  const double epsilon = std::numeric_limits<double>::epsilon();
  double diff = .01;

  auto val_a = [](int i, int j) { return 10 * i + j; };

  dca::linalg::Matrix<double, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);

  for (int sg : {1, -1})
    for (int ia : {0, 1, 4})
      for (int ja : {0, 2, 3}) {
        dca::linalg::Matrix<double, dca::linalg::CPU> b(a);
        b(ia, ja) += sg * diff;
        double err = std::abs(epsilon * b(ia, ja));

        EXPECT_NEAR(diff, dca::linalg::matrixop::difference(a, b, 2 * diff), err);
        EXPECT_NEAR(diff, dca::linalg::matrixop::difference(a, b, diff + err), err);
        auto diffcalc = dca::linalg::matrixop::difference(a, b, 2 * diff);
        EXPECT_NEAR(diff, dca::linalg::matrixop::difference(a, b, diffcalc), err);
        EXPECT_THROW(dca::linalg::matrixop::difference(a, b, diffcalc - err), std::logic_error);
      }
}

TEST(MatrixCPUTest, DeterminantAndLogDeterminant) {
  dca::linalg::Matrix<double, dca::linalg::CPU> m(3, 3);
  m(0, 0) = 3, m(0, 1) = -1, m(0, 2) = 0.5;
  m(1, 0) = 2;
  m(1, 1) = 2, m(1, 2) = -0.5;
  m(2, 0) = 4;
  m(2, 1) = -1, m(2, 2) = 0;

  const double det = dca::linalg::matrixop::determinant(m);
  const auto [log_det, sign] = dca::linalg::matrixop::logDeterminant(m);
  EXPECT_NEAR(-4.5, det, 1e-14);
  EXPECT_NEAR(std::log(std::abs(det)), log_det, 1e-14);
  EXPECT_EQ(det >= 0 ? 1 : -1, sign);

  dca::linalg::Matrix<double, dca::linalg::CPU> m2(2, 2);
  m2(0, 0) = 3, m2(0, 1) = 6;
  m2(1, 0) = 1;
  m2(1, 1) = 2;
  EXPECT_NEAR(0, dca::linalg::matrixop::determinantIP(m2), 1e-14);
}

template <typename ScalarType>
class MatrixopMixedCPUTest : public ::testing::Test {
public:
  static const ScalarType epsilon;
};
template <typename ScalarType>
const ScalarType MatrixopMixedCPUTest<ScalarType>::epsilon = std::numeric_limits<ScalarType>::epsilon();

TYPED_TEST_CASE(MatrixopMixedCPUTest, FloatingPointTypes);

TYPED_TEST(MatrixopMixedCPUTest, Gemm) {
  using ScalarType = TypeParam;
  using ComplexType = std::complex<ScalarType>;
  auto val_real = [](int i, int j) { return 3 * i - 2 * j; };
  auto val_complex = [](int i, int j) { return ComplexType(i - 2 * j * j, 4 * i - 3 * j); };
  auto val_c = [](int i, int j) { return ComplexType(2 * i - j, 2 * j - i); };

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

      if (transa != 'C') {
        dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
        dca::linalg::Matrix<ComplexType, dca::linalg::CPU> b(size_b);
        dca::linalg::Matrix<ComplexType, dca::linalg::CPU> c(size_c);

        testing::setMatrixElements(a, val_real);
        testing::setMatrixElements(b, val_complex);
        testing::setMatrixElements(c, val_c);

        dca::linalg::Matrix<ComplexType, dca::linalg::CPU> cref(c);

        testing::refGemm(transa, transb, ComplexType(1), a, b, ComplexType(0), cref);
        dca::linalg::matrixop::gemm(transa, transb, a, b, c);

        for (int j = 0; j < c.nrCols(); ++j)
          for (int i = 0; i < c.nrRows(); ++i) {
            EXPECT_GE(500 * this->epsilon, std::abs(cref(i, j) - c(i, j)));
          }
      }
      if (transb != 'C') {
        dca::linalg::Matrix<ComplexType, dca::linalg::CPU> a(size_a);
        dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size_b);
        dca::linalg::Matrix<ComplexType, dca::linalg::CPU> c(size_c);

        testing::setMatrixElements(a, val_complex);
        testing::setMatrixElements(b, val_real);
        testing::setMatrixElements(c, val_c);

        dca::linalg::Matrix<ComplexType, dca::linalg::CPU> cref(c);

        testing::refGemm(transa, transb, ComplexType(1), a, b, ComplexType(0), cref);
        dca::linalg::matrixop::gemm(transa, transb, a, b, c);

        for (int j = 0; j < c.nrCols(); ++j)
          for (int i = 0; i < c.nrRows(); ++i) {
            EXPECT_GE(500 * this->epsilon, std::abs(cref(i, j) - c(i, j)));
          }
      }
    }
  }
}

TYPED_TEST(MatrixopRealCPUTest, InverseAndDeterminant) {
  using ScalarType = TypeParam;
  dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 42);
  auto values = [&rng](int, int) { return rng(); };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(7);
  testing::setMatrixElements(mat, values);
  auto mat_copy(mat);
  auto mat_copy2(mat);

  dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
  dca::linalg::Vector<ScalarType, dca::linalg::CPU> work;

  const ScalarType det = dca::linalg::matrixop::inverseAndDeterminant(mat);
  dca::linalg::matrixop::inverse(mat_copy);

  const auto [log_det, sign] = dca::linalg::matrixop::inverseAndLogDeterminant(mat_copy2);

  const auto tolerance = 500 * this->epsilon;
  EXPECT_NEAR(dca::linalg::matrixop::determinant(mat_copy), det, tolerance);
  EXPECT_NEAR(std::log(std::abs(det)), log_det, tolerance);
  EXPECT_EQ(det >= 0 ? 1 : -1, sign);

  EXPECT_TRUE(dca::linalg::matrixop::areNear(mat_copy, mat, tolerance));
  EXPECT_TRUE(dca::linalg::matrixop::areNear(mat_copy2, mat, tolerance));
}
