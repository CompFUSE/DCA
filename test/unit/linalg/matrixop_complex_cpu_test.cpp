// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file tests the matrix operations with the Matrix<CPU> class with complex element type.

#include "dca/linalg/matrixop.hpp"

#include <array>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "gtest/gtest.h"

#include "dca/linalg/lapack/use_device.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/blas/blas3.hpp"
#include "cpu_test_util.hpp"
#include "matrixop_reference.hpp"

template <typename ComplexType>
class MatrixopComplexCPUTest : public ::testing::Test {
public:
  static const typename ComplexType::value_type epsilon;
};
template <typename ComplexType>
const typename ComplexType::value_type MatrixopComplexCPUTest<ComplexType>::epsilon =
    std::numeric_limits<typename ComplexType::value_type>::epsilon();

typedef ::testing::Types<std::complex<float>, std::complex<double>> ComplexFloatingPointTypes;
TYPED_TEST_CASE(MatrixopComplexCPUTest, ComplexFloatingPointTypes);

TYPED_TEST(MatrixopComplexCPUTest, Gemv) {
  using ScalarType = TypeParam;
  auto val_a = [](int i, int j) { return ScalarType(3 * i - 2 * j, 1 - i * i + j); };
  auto val_x = [](int i) { return ScalarType(i, 4 * i - 3 * i * i); };
  auto val_y = [](int i) { return ScalarType(2 * i, 2 * i * i); };

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
        EXPECT_GE(500 * this->epsilon, std::abs(yref[i] - y[i]));
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
        EXPECT_GE(500 * this->epsilon, std::abs(yref[i] - y[i]));
      }
    }
  }
}

TYPED_TEST(MatrixopComplexCPUTest, Gemm) {
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

      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> cref(c);

      testing::refGemm('N', 'N', ScalarType(1), a, b, ScalarType(0), cref);
      dca::linalg::matrixop::gemm(a, b, c);

      for (int j = 0; j < c.nrCols(); ++j)
        for (int i = 0; i < c.nrRows(); ++i) {
          EXPECT_GE(500 * this->epsilon, std::abs(cref(i, j) - c(i, j)));
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
          EXPECT_GE(500 * this->epsilon, std::abs(cref(i, j) - c(i, j)));
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
              EXPECT_GE(500 * this->epsilon, std::abs(cref(i, j) - c(i, j)));
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
              EXPECT_GE(500 * this->epsilon, std::abs(cref(i, j) - c(i, j)));
            }
        }
      }
    }
  }
}

TYPED_TEST(MatrixopComplexCPUTest, MultiplyDiagonal) {
  using ScalarType = TypeParam;
  std::pair<int, int> size_a(3, 5);
  auto val_a = [](int i, int j) { return ScalarType(3 * i - 2 * j, 1 - i * i + j); };
  auto val_d = [](int i) { return ScalarType(i, 1 - i); };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size_a);
  testing::setMatrixElements(a, val_a);
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b(size_a);
  {
    dca::linalg::Vector<ScalarType, dca::linalg::CPU> d(a.nrRows());
    testing::setVectorElements(d, val_d);

    dca::linalg::matrixop::multiplyDiagonalLeft(d, a, b);

    for (int j = 0; j < a.nrCols(); ++j)
      for (int i = 0; i < a.nrRows(); ++i) {
        EXPECT_GE(10 * this->epsilon, std::abs(b(i, j) - d[i] * a(i, j)));
      }
  }
  {
    dca::linalg::Vector<ScalarType, dca::linalg::CPU> d(a.nrCols());
    testing::setVectorElements(d, val_d);

    dca::linalg::matrixop::multiplyDiagonalRight(a, d, b);

    for (int j = 0; j < a.nrCols(); ++j)
      for (int i = 0; i < a.nrRows(); ++i) {
        EXPECT_GE(10 * this->epsilon, std::abs(b(i, j) - d[j] * a(i, j)));
      }
  }
}

TYPED_TEST(MatrixopComplexCPUTest, Trsm) {
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

      dca::linalg::Matrix<ScalarType, dca::linalg::CPU> b_orig(b);

      dca::linalg::matrixop::trsm(uplo, diag, a, b);
      dca::linalg::blas::trmm("L", &uplo, "N", &diag, m, n, ScalarType(1), a.ptr(),
                              a.leadingDimension(), b.ptr(), b.leadingDimension());

      for (int j = 0; j < b.nrCols(); ++j)
        for (int i = 0; i < b.nrRows(); ++i) {
          EXPECT_GE(500 * this->epsilon, std::abs(b_orig(i, j) - b(i, j)));
        }
    }
  }
}

TYPED_TEST(MatrixopComplexCPUTest, Eigensolver) {
  using ScalarType = TypeParam;
  int size = 2;
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size);
  mat(0, 0) = ScalarType(0, 2);
  mat(0, 1) = ScalarType(0, 0);
  mat(1, 0) = ScalarType(0, 1);
  mat(1, 1) = ScalarType(0, 1);

  dca::linalg::Vector<ScalarType, dca::linalg::CPU> w;
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> vl;
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> vr;

  dca::linalg::matrixop::eigensolver('V', 'V', mat, w, vl, vr);
  EXPECT_EQ(w.size(), size);
  EXPECT_EQ(vl.size(), std::make_pair(size, size));
  EXPECT_EQ(vr.size(), std::make_pair(size, size));

  EXPECT_GE(100 * this->epsilon, std::abs(ScalarType(0, 1) - w[0]));
  EXPECT_GE(100 * this->epsilon, std::abs(ScalarType(-1) - vl(0, 0) / vl(1, 0)));
  EXPECT_GE(100 * this->epsilon, std::abs(vr(0, 0)));
  EXPECT_GE(100 * this->epsilon, std::abs(ScalarType(1) - vr(1, 0)));

  EXPECT_GE(100 * this->epsilon, std::abs(ScalarType(0, 2) - w[1]));
  EXPECT_GE(100 * this->epsilon, std::abs(ScalarType(1) - vl(0, 1)));
  EXPECT_GE(100 * this->epsilon, std::abs(vl(1, 1)));
  EXPECT_GE(100 * this->epsilon, std::abs(ScalarType(1) - vr(0, 1) / vr(1, 1)));
}

TYPED_TEST(MatrixopComplexCPUTest, EigensolverHermitian) {
  using ScalarType = TypeParam;
  int size = 2;
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size);
  mat(0, 0) = ScalarType(2, 0);
  mat(0, 1) = ScalarType(0, 1);
  mat(1, 0) = ScalarType(0, -1);
  mat(1, 1) = ScalarType(2, 0);

  dca::linalg::Vector<typename ScalarType::value_type, dca::linalg::CPU> w;
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> v;

  dca::linalg::matrixop::eigensolverHermitian('V', 'U', mat, w, v);
  EXPECT_EQ(w.size(), size);
  EXPECT_EQ(v.size(), std::make_pair(size, size));

  EXPECT_NEAR(1., w[0], 100 * this->epsilon);
  EXPECT_GE(100 * this->epsilon, std::abs(ScalarType(0, -1) - v(0, 0) / v(1, 0)));

  EXPECT_NEAR(3., w[1], 100 * this->epsilon);
  EXPECT_GE(100 * this->epsilon, std::abs(ScalarType(0, 1) - v(0, 1) / v(1, 1)));
}

TYPED_TEST(MatrixopComplexCPUTest, Laset) {
  using ScalarType = TypeParam;
  std::pair<int, int> size(3, 35);
  ScalarType diag(3.4, 1.11);
  ScalarType offdiag(-1.4, 0.1);

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

TYPED_TEST(MatrixopComplexCPUTest, InsertRowCol) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2(3, 4);
  auto val = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };
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

      EXPECT_EQ(ScalarType(0), mat_test(ii, j));
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

      EXPECT_EQ(ScalarType(0), mat_test(i, jj));
    }
  }
}

TYPED_TEST(MatrixopComplexCPUTest, Inverse) {
  using ScalarType = TypeParam;
  int size = 6;
  auto val = [size](int i, int j) {
    return std::polar<typename ScalarType::value_type>(2. / ((4 + size + i - j) % size + 1), i + j);
  };
  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size);
  testing::setMatrixElements(mat, val);

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> invmat(mat);
  dca::linalg::matrixop::inverse(invmat);

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

TYPED_TEST(MatrixopComplexCPUTest, SmallInverse) {
  using ScalarType = TypeParam;

  for (int size = 1; size < 3; ++size) {
    auto val = [size](int i, int j) {
      return std::polar<typename ScalarType::value_type>(2. / ((4 + size + i - j) % size + 1), i + j);
    };
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> mat(size);
    testing::setMatrixElements(mat, val);

    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> invmat(mat);
    dca::linalg::matrixop::smallInverse(invmat);

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
}

TYPED_TEST(MatrixopComplexCPUTest, PseudoInverse) {
  using ScalarType = TypeParam;
  auto val = [](int i, int j) {
    return std::polar<typename ScalarType::value_type>(2. / ((5 + i - j) % 3 + 1), i + j);
  };
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
          EXPECT_GE(500 * this->epsilon, std::abs(ScalarType(1) - res(i, j)));
        else
          EXPECT_GE(500 * this->epsilon, std::abs(res(i, j)));
      }
    }

    // Check eigenvalue exclusion below threshold.
    testing::setMatrixElements(mat, val0);
    mat(0, 0) = ScalarType(20.);
    mat(1, 1) = ScalarType(.1);
    dca::linalg::matrixop::pseudoInverse(mat, invmat, 1e-6);
    if (size.first <= size.second)
      dca::linalg::matrixop::gemm(mat, invmat, res);
    else
      dca::linalg::matrixop::gemm(invmat, mat, res);

    for (int j = 0; j < res.nrCols(); ++j) {
      for (int i = 0; i < res.nrRows(); ++i) {
        if (i == j)
          EXPECT_GE(100 * this->epsilon, std::abs(ScalarType(1) - res(i, j)));
        else
          EXPECT_GE(100 * this->epsilon, std::abs(res(i, j)));
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
          EXPECT_GE(100 * this->epsilon, std::abs(ScalarType(1) - res(i, j)));
        else
          EXPECT_GE(100 * this->epsilon, std::abs(res(i, j)));
      }
    }
  }
}

TYPED_TEST(MatrixopComplexCPUTest, RemoveRowCol) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2(3, 4);
  auto val = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };
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

TYPED_TEST(MatrixopComplexCPUTest, CopyRow) {
  using ScalarType = TypeParam;
  const ScalarType checked(-1000);
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

  auto val_a = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };
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
        c(i_dests[i], j) = checked;
      }
      for (int i = 0; i < b.nrRows(); ++i)
        if (c(i, j) != checked)
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
        c(i_dests[i], j) = checked;
      }
      for (int i = 0; i < b.nrRows(); ++i)
        if (c(i, j) != checked)
          EXPECT_EQ(b(i, j), c(i, j));
    }
  }
}

TYPED_TEST(MatrixopComplexCPUTest, CopyCol) {
  using ScalarType = TypeParam;
  const ScalarType checked(-1000);
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

  auto val_a = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };
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
        c(i, j_dests[j]) = checked;
      }
      for (int j = 0; j < b.nrCols(); ++j)
        if (c(i, j) != checked)
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
        c(i, j_dests[j]) = checked;
      }
      for (int j = 0; j < b.nrCols(); ++j)
        if (c(i, j) != checked)
          EXPECT_EQ(b(i, j), c(i, j));
    }
  }
}

TYPED_TEST(MatrixopComplexCPUTest, ScaleRow) {
  using ScalarType = TypeParam;
  const ScalarType checked(-1000);
  std::pair<int, int> size2_a(4, 3);

  dca::linalg::Vector<int, dca::linalg::CPU> is(3);
  is[0] = 0;
  is[1] = 2;
  is[2] = 3;

  dca::linalg::Vector<ScalarType, dca::linalg::CPU> vals(3);
  vals[0] = 3.4;
  vals[1] = -1.2;
  vals[2] = 7.7;

  auto val_a = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };

  dca::linalg::Matrix<ScalarType, dca::linalg::CPU> a(size2_a);
  testing::setMatrixElements(a, val_a);

  {
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(a);

    for (int i = 0; i < is.size(); ++i)
      dca::linalg::matrixop::scaleRow(c, is[i], vals[i]);

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
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(a);

    dca::linalg::matrixop::scaleRows(c, is, vals);

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

TYPED_TEST(MatrixopComplexCPUTest, ScaleCol) {
  using ScalarType = TypeParam;
  const ScalarType checked(-1000);
  std::pair<int, int> size2_a(3, 4);

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
    dca::linalg::Matrix<ScalarType, dca::linalg::CPU> c(a);

    for (int j = 0; j < js.size(); ++j)
      dca::linalg::matrixop::scaleCol(c, js[j], vals[j]);

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

TYPED_TEST(MatrixopComplexCPUTest, SwapRow) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2_a(4, 3);

  auto val_a = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };

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

TYPED_TEST(MatrixopComplexCPUTest, SwapCol) {
  using ScalarType = TypeParam;
  std::pair<int, int> size2_a(3, 4);

  auto val_a = [](int i, int j) { return ScalarType(1 + i * i + j, 10 * i + j); };

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

TYPED_TEST(MatrixopComplexCPUTest, Multiply) {
  using ScalarType = typename TypeParam::value_type;

  auto val_a = [](int i, int j) { return std::complex<ScalarType>(1 + i + j, 10 * j * i); };
  auto val_b = [](int i, int j) { return std::complex<ScalarType>(i - j, j * j); };

  using RealMatrix = dca::linalg::Matrix<ScalarType, dca::linalg::CPU>;
  using CmplMatrix = dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU>;

  CmplMatrix a(std::make_pair(5, 5));
  CmplMatrix b(std::make_pair(5, 5));
  CmplMatrix c(std::make_pair(5, 5));

  testing::setMatrixElements(a, val_a);
  testing::setMatrixElements(b, val_b);

  std::array<RealMatrix, 2> a_split{RealMatrix(a.size()), RealMatrix(a.size())};
  std::array<RealMatrix, 2> b_split{RealMatrix(b.size()), RealMatrix(b.size())};
  std::array<RealMatrix, 2> c_split{RealMatrix(c.size()), RealMatrix(c.size())};

  auto split_real_imag = [](const CmplMatrix& m, std::array<RealMatrix, 2>& m_split) {
    for (int j = 0; j < m.nrCols(); ++j)
      for (int i = 0; i < m.nrRows(); ++i) {
        m_split[0](i, j) = m(i, j).real();
        m_split[1](i, j) = m(i, j).imag();
      }
  };
  split_real_imag(a, a_split);
  split_real_imag(b, b_split);

  std::array<RealMatrix, 5> work_space;
  for (char transa : {'N', 'T', 'C'})
    for (char transb : {'N', 'T', 'C'}) {
      // Compute result with direct multiplication.
      dca::linalg::matrixop::gemm(transa, transb, a, b, c);
      // compute with 3M algorithm,
      dca::linalg::matrixop::multiply(transa, transb, a_split, b_split, c_split, work_space);
      // Confront.
      const ScalarType tolerance = 500 * this->epsilon;
      for (int j = 0; j < c.nrCols(); ++j)
        for (int i = 0; i < c.nrRows(); ++i) {
          EXPECT_NEAR(c(i, j).real(), c_split[0](i, j), tolerance);
          EXPECT_NEAR(c(i, j).imag(), c_split[1](i, j), tolerance);
        }
    }
}

TYPED_TEST(MatrixopComplexCPUTest, MultiplyRealArg) {
  using ScalarType = typename TypeParam::value_type;

  auto val_a = [](int i, int j) { return std::complex<ScalarType>(i + j * i, -0.6 * j); };
  auto val_b = [](int i, int j) { return std::complex<ScalarType>(i - 2. * j, 0.); };

  using RealMatrix = dca::linalg::Matrix<ScalarType, dca::linalg::CPU>;
  using CmplMatrix = dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU>;

  CmplMatrix a(std::make_pair(7, 7));
  CmplMatrix b(std::make_pair(7, 7));
  CmplMatrix c(std::make_pair(7, 7));

  testing::setMatrixElements(a, val_a);
  testing::setMatrixElements(b, val_b);

  std::array<RealMatrix, 2> a_split{RealMatrix(a.size()), RealMatrix(a.size())};
  RealMatrix b_re(b.size());
  std::array<RealMatrix, 2> c_split{RealMatrix(c.size()), RealMatrix(c.size())};

  for (int j = 0; j < a.nrCols(); ++j)
    for (int i = 0; i < a.nrRows(); ++i) {
      a_split[0](i, j) = a(i, j).real();
      a_split[1](i, j) = a(i, j).imag();
    }
  for (int j = 0; j < b.nrCols(); ++j)
    for (int i = 0; i < b.nrRows(); ++i)
      b_re(i, j) = b(i, j).real();

  for (char transa : {'N', 'T'})
    for (char transb : {'N', 'T'}) {
      // Compute result with direct multiplication.
      dca::linalg::matrixop::gemm(transa, transb, a, b, c);
      // compute with 3M algorithm,
      dca::linalg::matrixop::multiply(transa, transb, a_split, b_re, c_split);
      // Confront.
      const ScalarType tolerance = 500 * this->epsilon;
      for (int j = 0; j < c.nrCols(); ++j)
        for (int i = 0; i < c.nrRows(); ++i) {
          EXPECT_NEAR(c(i, j).real(), c_split[0](i, j), tolerance);
          EXPECT_NEAR(c(i, j).imag(), c_split[1](i, j), tolerance);
        }
    }
}
