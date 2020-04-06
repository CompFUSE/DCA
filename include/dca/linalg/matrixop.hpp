// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides the matrix interface for the following matrix operations:
// - copyCol, copyRow, copyCols, copyRows
// - difference
// - real
// - insertCol, insertRow (for CPU matrices only)
// - inverse
// - inverseAndDeterminant
// - removeCol, removeCols, removeRow, removeRows, removeRowAndCol, removeRowAndCols
// - scaleCol, scaleRow, scaleRows
// - swapCol, swapRow, swapRowAndCol
// - swapCols, swapRows (for GPU matrices only)
// - gemm
// - multiply
// - trsm
// - determinant
// - logDeterminant
// - eigensolver (non-symmetric / symmetric / Hermitian)
// - pseudoInverse

#ifndef DCA_LINALG_MATRIXOP_HPP
#define DCA_LINALG_MATRIXOP_HPP

#include <cassert>
#include <cstring>
#include <tuple>

#include "dca/linalg/blas/use_device.hpp"
#include "dca/linalg/lapack/use_device.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/util/util_lapack.hpp"
#include "dca/linalg/util/util_matrixop.hpp"
#include "dca/linalg/vector.hpp"

#ifdef DCA_HAVE_CUDA
#include "dca/linalg/blas/kernels_gpu.hpp"
#endif

namespace dca {
namespace linalg {
namespace matrixop {
// dca::linalg::matrixop::

// Copies the matrix mat in a.
// Preconditions: lda >= mat.nrRows().
template <typename Scalar>
inline void copyMatrixToArray(const Matrix<Scalar, CPU>& mat, Scalar* a, int lda) {
  assert(lda >= mat.nrRows());
  lapack::lacpy("A", mat.nrRows(), mat.nrCols(), mat.ptr(), mat.leadingDimension(), a, lda);
}

// Copies the m by n matrix stored in a to the matrix mat.
// Preconditions: lda >= m.
template <typename Scalar>
inline void copyArrayToMatrix(int m, int n, const Scalar* a, int lda, Matrix<Scalar, CPU>& mat) {
  assert(lda >= m);
  mat.resizeNoCopy(std::make_pair(m, n));
  lapack::lacpy("A", mat.nrRows(), mat.nrCols(), a, lda, mat.ptr(), mat.leadingDimension());
}

// Copies the jx-th column of mat_x into the jy-th column of mat_y.
// In/Out: mat_y
// Preconditions: mat_x.nrRows() == mat_y.nrRows(),
//                0 <= jx < mat_x.nrCols(), 0 <= jy < mat_y.nrCols().
template <typename Scalar, DeviceType device_name>
inline void copyCol(const Matrix<Scalar, device_name>& mat_x, int jx,
                    Matrix<Scalar, device_name>& mat_y, int jy, int thread_id = 0, int stream_id = 0) {
  assert(jx >= 0 && jx < mat_x.nrCols());
  assert(jy >= 0 && jy < mat_y.nrCols());
  assert(mat_x.nrRows() == mat_y.nrRows());

  blas::UseDevice<device_name>::copy(mat_x.nrRows(), mat_x.ptr(0, jx), 1, mat_y.ptr(0, jy), 1,
                                     thread_id, stream_id);
}

// Copies the j_x[i]-th column of mat_x into the j_y[i]-th column of mat_y, for 0 <= i < j_x.size().
// In/Out: mat_y
// Preconditions: j_x.size() <= j_y.size(), mat_x.nrRows() == mat_y.nrRows()
//                0 <= j_x[i] < mat_x.nrCols() for 0 <= i < j_x.size(),
//                0 <= j_y[i] < mat_y.nrCols() for 0 <= i < j_x.size().
template <typename Scalar, class Vec>
inline void copyCols(const Matrix<Scalar, CPU>& mat_x, const Vec& j_x, Matrix<Scalar, CPU>& mat_y,
                     const Vec& j_y, int /*thread_id*/ = 0, int /*stream_id*/ = 0) {
  assert(j_x.size() <= j_y.size());

  for (int ind_j = 0; ind_j < j_x.size(); ++ind_j)
    copyCol(mat_x, j_x[ind_j], mat_y, j_y[ind_j]);
}
#ifdef DCA_HAVE_CUDA
template <typename Scalar>
inline void copyCols(const Matrix<Scalar, GPU>& mat_x, const Vector<int, GPU>& j_x,
                     Matrix<Scalar, GPU>& mat_y, const Vector<int, GPU>& j_y, int thread_id = 0,
                     int stream_id = 0) {
  assert(j_x.size() <= j_y.size());
  assert(mat_x.nrRows() == mat_y.nrRows());

  blas::copyCols(mat_x.nrRows(), j_x.size(), j_x.ptr(), mat_x.ptr(), mat_x.leadingDimension(),
                 j_y.ptr(), mat_y.ptr(), mat_y.leadingDimension(), thread_id, stream_id);
}

// Copies the j_x columns of mat_x into the  mat_y, for 0 <= i < j_x.size().
// In/Out: mat_y
// Preconditions: mat_x.nrRows() == mat_y.nrRows()
//                0 <= j_x[i] < mat_x.nrCols() for 0 <= i < j_x.size(),
template <typename Scalar>
inline void copyCols(const Matrix<Scalar, GPU>& mat_x, const Vector<int, GPU>& j_x,
                     Matrix<Scalar, GPU>& mat_y, int thread_id = 0, int stream_id = 0) {
  assert(mat_x.nrRows() == mat_y.nrRows());

  blas::copyCols(mat_x.nrRows(), j_x.size(), j_x.ptr(), mat_x.ptr(), mat_x.leadingDimension(),
                 mat_y.ptr(), mat_y.leadingDimension(), thread_id, stream_id);
}
#endif  // DCA_HAVE_CUDA

// Copies the ix-th row of mat_x into the iy-th row of mat_y.
// In/Out: mat_y
// Preconditions: mat_x.nrCols() == mat_y.nrCols(),
//                0 <= ix < mat_x.nrRows(), 0 <= iy < mat_y.nrRows().
template <typename Scalar, DeviceType device_name>
inline void copyRow(const Matrix<Scalar, device_name>& mat_x, int ix,
                    Matrix<Scalar, device_name>& mat_y, int iy, int thread_id = 0, int stream_id = 0) {
  assert(ix >= 0 && ix < mat_x.nrRows());
  assert(iy >= 0 && iy < mat_y.nrRows());
  assert(mat_x.nrCols() == mat_y.nrCols());

  blas::UseDevice<device_name>::copy(mat_x.nrCols(), mat_x.ptr(ix, 0), mat_x.leadingDimension(),
                                     mat_y.ptr(iy, 0), mat_y.leadingDimension(), thread_id,
                                     stream_id);
}

// Copies the i_x[i]-th row of mat_x into the i_y[i]-th row of mat_y, for 0 <= i < i_x.size().
// In/Out: mat_y
// Preconditions: i_x.size() <= i_y.size(), mat_x.nrCols() == mat_y.nrCols()
//                0 <= i_x[i] < mat_x.nrRows() for 0 <= i < i_x.size(),
//                0 <= i_y[i] < mat_y.nrRows() for 0 <= i < i_x.size().
template <typename Scalar, class Vec>
inline void copyRows(const Matrix<Scalar, CPU>& mat_x, const Vec& i_x, Matrix<Scalar, CPU>& mat_y,
                     const Vec& i_y, int /*thread_id*/ = 0, int /*stream_id*/ = 0) {
  assert(i_x.size() <= i_y.size());
  assert(mat_x.nrCols() == mat_y.nrCols());

  for (int j = 0; j < mat_x.nrCols(); ++j)
    for (int ind_i = 0; ind_i < i_x.size(); ++ind_i)
      mat_y(i_y[ind_i], j) = mat_x(i_x[ind_i], j);
}
#ifdef DCA_HAVE_CUDA
template <typename Scalar>
inline void copyRows(const Matrix<Scalar, GPU>& mat_x, const Vector<int, GPU>& i_x,
                     Matrix<Scalar, GPU>& mat_y, const Vector<int, GPU>& i_y, int thread_id = 0,
                     int stream_id = 0) {
  assert(i_x.size() <= i_y.size());
  assert(mat_x.nrCols() == mat_y.nrCols());

  blas::copyRows(mat_x.nrCols(), i_x.size(), i_x.ptr(), mat_x.ptr(), mat_x.leadingDimension(),
                 i_y.ptr(), mat_y.ptr(), mat_y.leadingDimension(), thread_id, stream_id);
}

// Copies the i_x rows of mat_x into  mat_y, for 0 <= i < i_x.size().
// In/Out: mat_y
// Preconditions: mat_x.nrCols() == mat_y.nrCols()
//                0 <= i_x[i] < mat_x.nrRows() for 0 <= i < i_x.size().
template <typename Scalar>
inline void copyRows(const Matrix<Scalar, GPU>& mat_x, const Vector<int, GPU>& i_x,
                     Matrix<Scalar, GPU>& mat_y, int thread_id = 0, int stream_id = 0) {
  assert(mat_x.nrCols() == mat_y.nrCols());

  blas::copyRows(mat_x.nrCols(), i_x.size(), i_x.ptr(), mat_x.ptr(), mat_x.leadingDimension(),
                 mat_y.ptr(), mat_y.leadingDimension(), thread_id, stream_id);
}
#endif  // DCA_HAVE_CUDA

// Returns the difference of two matrices in terms of max_i,j(|a(i, j) - b(i, j)|).
// If the difference is larger than the threshold a std::logic_error exception is thrown,
// and if NDEBUG is not defined each difference which exceeds the threshold is printed.
// Preconditions: a.size() == b.size().
template <typename Scalar>
auto difference(const Matrix<Scalar, CPU>& a, const Matrix<Scalar, CPU>& b,
                double diff_threshold = 1e-3) {
  assert(a.size() == b.size());

  auto max_diff = std::abs(Scalar(0));

  for (int j = 0; j < a.nrCols(); ++j) {
    for (int i = 0; i < a.nrRows(); ++i) {
      max_diff = std::max(max_diff, std::abs(a(i, j) - b(i, j)));
    }
  }

  if (max_diff > diff_threshold) {
#ifndef NDEBUG
    std::stringstream s;
    for (int i = 0; i < a.nrRows(); ++i) {
      for (int j = 0; j < a.nrCols(); ++j) {
        if (std::abs(a(i, j) - b(i, j)) <= diff_threshold)
          s << 0. << "\t";
        else
          s << a(i, j) - b(i, j) << "\t";
      }
      s << "\n";
    }
    s << std::endl;
    std::cout << s.str();
#endif  // NDEBUG

    throw std::logic_error(__FUNCTION__);
  }

  return max_diff;
}
template <typename Scalar, DeviceType device_name>
auto difference(const Matrix<Scalar, device_name>& a, const Matrix<Scalar, CPU>& b,
                double diff_threshold = 1e-3) {
  Matrix<Scalar, CPU> cp_a(a);
  return difference(cp_a, b, diff_threshold);
}
template <typename Scalar, DeviceType device_name_a, DeviceType device_name_b>
auto difference(const Matrix<Scalar, device_name_a>& a, const Matrix<Scalar, device_name_b>& b,
                double diff_threshold = 1e-3) {
  Matrix<Scalar, CPU> cp_b(b);
  return difference(a, cp_b, diff_threshold);
}

// Returns the real part of a matrix.
// In: a
// TODO test.
template <typename Scalar>
Matrix<Scalar, CPU> real(const Matrix<std::complex<Scalar>, CPU>& a) {
  Matrix<Scalar, CPU> a_re(a.size());
  for (int j = 0; j < a.nrCols(); ++j)
    for (int i = 0; i < a.nrRows(); ++i)
      a_re(i, j) = std::real(a(i, j));
  return a_re;
}

// Insert a column at position j. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= j < mat.nrCols() + 1.
// Postconditions: The elements of the inserted column are set to 0.
template <typename Scalar>
void insertCol(Matrix<Scalar, CPU>& mat, int j) {
  assert(j >= 0 && j < mat.nrCols() + 1);

  mat.resize(std::make_pair(mat.nrRows(), mat.nrCols() + 1));

  if (mat.nrRows() > 0 && j < mat.nrCols() - 1)
    memmove(mat.ptr(0, j + 1), mat.ptr(0, j),
            sizeof(Scalar) * (mat.nrCols() - 1 - j) * mat.leadingDimension());

  for (int i = 0; i < mat.nrRows(); ++i)
    mat(i, j) = 0;
}

// Insert a row at position i. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= i < mat.nrRows() + 1.
// Postconditions: The elements of the inserted row are set to 0.
template <typename Scalar>
void insertRow(Matrix<Scalar, CPU>& mat, int i) {
  assert(i >= 0 && i < mat.nrRows() + 1);

  mat.resize(std::make_pair(mat.nrRows() + 1, mat.nrCols()));

  if (i < mat.nrRows() - 1)
    for (int j = 0; j < mat.nrCols(); ++j)
      memmove(mat.ptr(i + 1, j), mat.ptr(i, j), sizeof(Scalar) * (mat.nrRows() - 1 - i));

  for (int j = 0; j < mat.nrCols(); ++j)
    mat(i, j) = 0;
}

// Computes the inverse of the matrix using the LU factorization.
// In/Out: mat
// Out: ipiv, work
// Preconditions: mat is a square matrix.
// Postconditions: ipiv and work are resized to the needed dimension.
template <typename Scalar, DeviceType device_name, template <typename, DeviceType> class MatrixType>
void inverse(MatrixType<Scalar, device_name>& mat, Vector<int, CPU>& ipiv,
             Vector<Scalar, device_name>& work) {
  assert(mat.is_square());

  ipiv.resizeNoCopy(mat.nrRows());

  lapack::UseDevice<device_name>::getrf(mat.nrRows(), mat.nrCols(), mat.ptr(),
                                        mat.leadingDimension(), ipiv.ptr());
  // Get optimal worksize.
  int lwork = util::getInverseWorkSize(mat);
  work.resizeNoCopy(lwork);

  lapack::UseDevice<device_name>::getri(mat.nrRows(), mat.ptr(), mat.leadingDimension(), ipiv.ptr(),
                                        work.ptr(), lwork);
}
template <typename Scalar, DeviceType device_name, template <typename, DeviceType> class MatrixType>
void inverse(MatrixType<Scalar, device_name>& mat) {
  Vector<int, CPU> ipiv;
  Vector<Scalar, device_name> work;
  inverse(mat, ipiv, work);
}

template <typename Scalar>
void smallInverse(Matrix<Scalar, CPU>& m_inv, Vector<int, CPU>& ipiv, Vector<Scalar, CPU>& work) {
  assert(m_inv.is_square());
  switch (m_inv.nrCols()) {
    case 1:
      m_inv(0, 0) = Scalar(1.) / m_inv(0, 0);
      break;
    case 2: {
      const Scalar det = m_inv(0, 0) * m_inv(1, 1) - m_inv(0, 1) * m_inv(1, 0);

      std::swap(m_inv(0, 0), m_inv(1, 1));
      m_inv(0, 0) /= det;
      m_inv(1, 1) /= det;
      std::swap(m_inv(1, 0), m_inv(0, 1));
      m_inv(1, 0) /= -det;
      m_inv(0, 1) /= -det;
      break;
    }
    case 3: {
      const Matrix<Scalar, CPU> m(m_inv);
      const Scalar det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
                         m(1, 0) * (m(0, 1) * m(2, 2) - m(0, 2) * m(2, 1)) +
                         m(2, 0) * (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1));
      m_inv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) / det;
      m_inv(0, 1) = -(m(0, 1) * m(2, 2) - m(0, 2) * m(2, 1)) / det;
      m_inv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) / det;
      m_inv(1, 0) = -(m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) / det;
      m_inv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) / det;
      m_inv(1, 2) = -(m(0, 0) * m(1, 2) - m(1, 0) * m(0, 2)) / det;
      m_inv(2, 0) = (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0)) / det;
      m_inv(2, 1) = -(m(0, 0) * m(2, 1) - m(0, 1) * m(2, 0)) / det;
      m_inv(2, 2) = (m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0)) / det;
      break;
    }
    default:
      inverse(m_inv, ipiv, work);
  }
}

template <class Scalar>
void smallInverse(Matrix<Scalar, CPU>& m_inv) {
  Vector<int, CPU> ipiv;
  Vector<Scalar, CPU> work;
  smallInverse(m_inv, ipiv, work);
}

// Computes in place the inverse of mat and the determinant of the inverse.
// In/Out: mat
// Returns: the determinant of mat^-1
// Precondition: mat is a non-singular real matrix.
template <typename Scalar, template <typename, DeviceType> class MatrixType>
Scalar inverseAndDeterminant(MatrixType<Scalar, CPU>& mat) {
  assert(mat.is_square());
  std::vector<int> ipiv(mat.nrRows());

  lapack::UseDevice<CPU>::getrf(mat.nrRows(), mat.nrCols(), mat.ptr(), mat.leadingDimension(),
                                ipiv.data());

  Scalar det = 1;
  for (int i = 0; i < mat.nrCols(); ++i) {
    det *= mat(i, i);
    if (ipiv[i] != i + 1)
      det *= -1;
  }

  const int lwork = util::getInverseWorkSize(mat);
  std::vector<Scalar> work(lwork);
  lapack::UseDevice<CPU>::getri(mat.nrRows(), mat.ptr(), mat.leadingDimension(), ipiv.data(),
                                work.data(), lwork);

  return 1. / det;
}

// Remove the j-th column. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= j < mat.nrCols().
template <typename Scalar>
void removeCol(Matrix<Scalar, CPU>& mat, int j) {
  assert(j >= 0 && j < mat.nrCols());

  if (mat.nrRows() > 0 && j < mat.nrCols() - 1)
    memmove(mat.ptr(0, j), mat.ptr(0, j + 1),
            sizeof(Scalar) * (mat.nrCols() - j - 1) * mat.leadingDimension());

  mat.resize(std::make_pair(mat.nrRows(), mat.nrCols() - 1));
}

#ifdef DCA_HAVE_CUDA
template <typename Scalar>
void removeCol(Matrix<Scalar, GPU>& mat, int j) {
  assert(j >= 0 && j < mat.nrCols());

  if (mat.nrRows() > 0 && j < mat.nrCols() - 1)
    blas::moveLeft(mat.nrRows(), mat.nrCols() - j, mat.ptr(0, j), mat.leadingDimension());

  mat.resize(std::make_pair(mat.nrRows(), mat.nrCols() - 1));
}
#endif  // DCA_HAVE_CUDA

// Remove columns in range [first, last]. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= first, last < mat.nrCols().
template <typename Scalar>
void removeCols(Matrix<Scalar, CPU>& mat, int first, int last) {
  const int n_removed = last - first + 1;
  const int n = mat.nrRows();
  const int m = mat.nrCols();
  assert(last < m and last >= first and first >= 0);

  if (n > 0 and last < m - 1)
    std::memmove(mat.ptr(0, first), mat.ptr(0, last + 1),
                 mat.leadingDimension() * (m - last - 1) * sizeof(Scalar));

  mat.resize(std::make_pair(n, m - n_removed));
}

// Remove the i-th row. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= i < mat.nrRows().
template <typename Scalar>
void removeRow(Matrix<Scalar, CPU>& mat, int i) {
  assert(i >= 0 && i < mat.nrRows());

  if (i < mat.nrRows() - 1)
    for (int j = 0; j < mat.nrCols(); ++j)
      memmove(mat.ptr(i, j), mat.ptr(i + 1, j), sizeof(Scalar) * (mat.nrRows() - i - 1));

  mat.resize(std::make_pair(mat.nrRows() - 1, mat.nrCols()));
}

#ifdef DCA_HAVE_CUDA
template <typename Scalar>
void removeRow(Matrix<Scalar, GPU>& mat, int i) {
  assert(i >= 0 && i < mat.nrRows());

  if (mat.nrCols() > 0 && i < mat.nrRows() - 1)
    blas::moveUp(mat.nrRows() - i, mat.nrCols(), mat.ptr(i, 0), mat.leadingDimension());

  mat.resize(std::make_pair(mat.nrRows() - 1, mat.nrCols()));
}
#endif  // DCA_HAVE_CUDA

// Remove rows in range [first, last]. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= first, last < mat.nrRows().
template <typename Scalar>
void removeRows(Matrix<Scalar, CPU>& mat, int first, int last) {
  const int n_removed = last - first + 1;
  const int n = mat.nrRows();
  const int m = mat.nrCols();
  assert(last < n and last >= first and first >= 0);

  if (last < n - 1)
    for (int j = 0; j < m; ++j)
      std::memmove(mat.ptr(first, j), mat.ptr(last + 1, j), (n - last - 1) * sizeof(Scalar));

  mat.resize(std::make_pair(n - n_removed, m));
}

// Remove the i-th row and the j-th column. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= i < mat.nrRows(), 0 <= j < mat.nrCols().
template <typename Scalar, DeviceType device_name>
inline void removeRowAndCol(Matrix<Scalar, device_name>& mat, int i, int j) {
  removeRow(mat, i);
  removeCol(mat, j);
}

// Remove the i-th row and the i-th column. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= i < mat.nrRows(), i < mat.nrCols().
template <typename Scalar, DeviceType device_name>
inline void removeRowAndCol(Matrix<Scalar, device_name>& mat, int i) {
  removeRowAndCol(mat, i, i);
}

// Remove rows and columns in range [first, last]. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= first, last < min(mat.nrRows(), mat.nrCols()).
template <typename Scalar>
void removeRowsAndCols(Matrix<Scalar, CPU>& mat, int first, int last) {
  removeCols(mat, first, last);
  removeRows(mat, first, last);
}

// Scales the j-th column of mat by val.
// In/Out: mat
// Preconditions: 0 <= j < mat.nrCols().
template <typename Scalar, DeviceType device_name>
inline void scaleCol(Matrix<Scalar, device_name>& mat, int j, Scalar val, int thread_id = 0,
                     int stream_id = 0) {
  assert(j >= 0 && j < mat.nrCols());
  blas::UseDevice<device_name>::scal(mat.nrRows(), val, mat.ptr(0, j), 1, thread_id, stream_id);
}

// Scales the i-th row of mat by val.
// In/Out: mat
// Preconditions: 0 <= i < mat.nrRow().
template <typename Scalar, DeviceType device_name>
inline void scaleRow(Matrix<Scalar, device_name>& mat, int i, Scalar val, int thread_id = 0,
                     int stream_id = 0) {
  assert(i >= 0 && i < mat.nrRows());
  blas::UseDevice<device_name>::scal(mat.nrCols(), val, mat.ptr(i, 0), mat.leadingDimension(),
                                     thread_id, stream_id);
}

// Scales the i[k]-th row of mat by val[k] for 0 <= k < i.size().
// In/Out: mat
// Preconditions: i.size() == val.size(), 0 <= i[k] < mat.nrRow() for 0 <= k < i.size().
template <typename Scalar>
inline void scaleRows(Matrix<Scalar, CPU>& mat, const Vector<int, CPU>& i,
                      const Vector<Scalar, CPU>& val, int /*thread_id*/ = 0, int /*stream_id*/ = 0) {
  assert(i.size() == val.size());

  for (int j = 0; j < mat.nrCols(); ++j)
    for (int ind = 0; ind < i.size(); ++ind)
      mat(i[ind], j) *= val[ind];
}
#ifdef DCA_HAVE_CUDA
template <typename Scalar>
inline void scaleRows(Matrix<Scalar, GPU>& mat, const Vector<int, GPU>& i,
                      const Vector<Scalar, GPU>& val, int thread_id = 0, int stream_id = 0) {
  assert(i.size() == val.size());

  blas::scaleRows(mat.nrCols(), i.size(), i.ptr(), val.ptr(), mat.ptr(), mat.leadingDimension(),
                  thread_id, stream_id);
}
#endif  // DCA_HAVE_CUDA

// Swaps the j1-th column with the j2-th column of mat.
// In/Out: mat
// Preconditions: 0 <= j1 < mat.nrCols(), 0 <= j2 < mat_y.nrCols().
template <typename Scalar, DeviceType device_name>
inline void swapCol(Matrix<Scalar, device_name>& mat, int j1, int j2, int thread_id = 0,
                    int stream_id = 0) {
  assert(j1 >= 0 && j1 < mat.nrCols());
  assert(j2 >= 0 && j2 < mat.nrCols());
  blas::UseDevice<device_name>::swap(mat.nrRows(), mat.ptr(0, j1), 1, mat.ptr(0, j2), 1, thread_id,
                                     stream_id);
}

// Swaps the j_1[i]-th column with the j_2[i]-th column of mat, for 0 <= i < j_1.size().
// In/Out: mat
// Preconditions: j_1.size() <= j_2.size()
//                0 <= j_1[i] < mat.nrCols() for 0 <= i < j_1.size(),
//                0 <= j_2[i] < mat.nrCols() for 0 <= i < j_1.size().
//                j_1[i] != j_1[j] for i != j, j_2[i] != j_2[j] for i != j,
//                j_1[i] != j_2[j] for all i, j.
#ifdef DCA_HAVE_CUDA
template <typename Scalar>
inline void swapCols(Matrix<Scalar, GPU>& mat, const Vector<int, GPU>& j_1,
                     const Vector<int, GPU>& j_2, int thread_id = 0, int stream_id = 0) {
  assert(j_1.size() <= j_2.size());
  blas::swapCols(mat.nrRows(), j_1.size(), j_1.ptr(), j_2.ptr(), mat.ptr(), mat.leadingDimension(),
                 thread_id, stream_id);
}
#endif  // DCA_HAVE_CUDA

// Swaps the i1-th row with the i2-th row of mat.
// In/Out: mat
// Preconditions: 0 <= i1 < mat.nrRows(), 0 <= i2 < mat_y.nrRows().
template <typename Scalar, DeviceType device_name>
inline void swapRow(Matrix<Scalar, device_name>& mat, int i1, int i2, int thread_id = 0,
                    int stream_id = 0) {
  assert(i1 >= 0 && i1 < mat.nrRows());
  assert(i2 >= 0 && i2 < mat.nrRows());
  blas::UseDevice<device_name>::swap(mat.nrCols(), mat.ptr(i1, 0), mat.leadingDimension(),
                                     mat.ptr(i2, 0), mat.leadingDimension(), thread_id, stream_id);
}

// Swaps the i_1[i]-th row with the i_2[i]-th row of mat, for 0 <= i < i_1.size().
// In/Out: mat
// Preconditions: i_2.size() == i_2.size()
//                0 <= i_1[i] < mat.nrRows() for 0 <= i < i_1.size(),
//                0 <= i_2[i] < mat.nrRows() for 0 <= i < i_1.size().
//                i_1[i] != i_1[j] for i != j, i_2[i] != i_2[j] for i != j,
//                i_1[i] != i_2[j] for all i, j.
#ifdef DCA_HAVE_CUDA
template <typename Scalar>
inline void swapRows(Matrix<Scalar, GPU>& mat, const Vector<int, GPU>& i_1,
                     const Vector<int, GPU>& i_2, int thread_id = 0, int stream_id = 0) {
  assert(i_1.size() == i_2.size());
  blas::swapRows(mat.nrCols(), i_1.size(), i_1.ptr(), i_2.ptr(), mat.ptr(), mat.leadingDimension(),
                 thread_id, stream_id);
}
#endif  // DCA_HAVE_CUDA

// Swaps the i1-th row with the i2-th row and the i1-th column with the i2-th column of mat.
// In/Out: mat
// Preconditions: 0 <= i1 < mat.nrRows(), i1 < mat.nrCols(),
//                0 <= i2 < mat.nrRows(), i2 < mat.nrCols().
template <typename Scalar, DeviceType device_name>
inline void swapRowAndCol(Matrix<Scalar, device_name>& mat, int i1, int i2, int thread_id = 0,
                          int stream_id = 0) {
  swapRow(mat, i1, i2, thread_id, stream_id);
  swapCol(mat, i1, i2, thread_id, stream_id);
}

// Performs the matrix-vector multiplication y <- alpha * op(a) * x + beta * y,
// where op(X) = X if transX == 'N', op(X) = transposed(X) if transX == 'T', and
// op(X) == conjugate_transposed(X) if transX == 'C' (X = a).
// In/Out: y ('In' only if beta != 0)
// Preconditions: transa should be one of the following: 'N', 'T' or 'C',
//                a.nrRows() == y.size() if transa == 'N', a.nrCols() == y.size() otherwise,
//                a.nrCols() == x.size() if transa == 'N', a.nrRows() == x.size() otherwise.
template <typename Scalar>
void gemv(char transa, Scalar alpha, const Matrix<Scalar, CPU>& a, const Vector<Scalar, CPU>& x,
          Scalar beta, Vector<Scalar, CPU>& y) {
  if (transa == 'N') {
    assert(a.nrRows() == y.size());
    assert(a.nrCols() == x.size());
  }
  else {
    assert(a.nrRows() == x.size());
    assert(a.nrCols() == y.size());
  }

  int lda = a.leadingDimension();

  blas::gemv(&transa, a.nrRows(), a.nrCols(), alpha, a.ptr(), lda, x.ptr(), 1, beta, y.ptr(), 1);
}

// Performs the matrix-vector multiplication y <- op(a) * x,
// where op(X) = X if transX == 'N', op(X) = transposed(X) if transX == 'T', and
// op(X) == conjugate_transposed(X) if transX == 'C' (X = a).
// Out: y
// Preconditions: transa should be one of the following: 'N', 'T' or 'C',
//                a.nrRows() == y.size() if transa == 'N', a.nrCols() == y.size() otherwise,
//                a.nrCols() == x.size() if transa == 'N', a.nrRows() == x.size() otherwise.
template <typename Scalar>
void gemv(char transa, const Matrix<Scalar, CPU>& a, const Vector<Scalar, CPU>& x,
          Vector<Scalar, CPU>& y) {
  gemv<Scalar>(transa, 1., a, x, 0., y);
}

// Performs the matrix-matrix multiplication c <- alpha * op(a) * op(b) + beta * c,
// where op(X) = X if transX == 'N', op(X) = transposed(X) if transX == 'T', and
// op(X) == conjugate_transposed(X) if transX == 'C' (X = a, b).
// In/Out: c ('In' only if beta != 0)
// Preconditions: transa and transb should be one of the following: 'N', 'T' or 'C',
//                a.nrRows() == c.nrRows() if transa == 'N', a.nrCols() == c.nrRows() otherwise,
//                b.nrCols() == c.nrCols() if transb == 'N', b.nrRows() == c.nrCols() otherwise,
//                ka == kb, where ka = a.nrCols() if transa == 'N', ka = a.nrRows() otherwise and
//                          kb = b.nrRows() if transb == 'N', kb = b.nrCols() otherwise.
template <typename Scalar, DeviceType device_name, template <typename, DeviceType> class MatrixA,
          template <typename, DeviceType> class MatrixB, template <typename, DeviceType> class MatrixC>
void gemm(char transa, char transb, Scalar alpha, const MatrixA<Scalar, device_name>& a,
          const MatrixB<Scalar, device_name>& b, Scalar beta, MatrixC<Scalar, device_name>& c,
          int thread_id = 0, int stream_id = 0) {
  int m = c.nrRows();
  int n = c.nrCols();
  int k;

  if (transa == 'N') {
    assert(a.nrRows() == m);
    k = a.nrCols();
  }
  else {
    assert(a.nrCols() == m);
    k = a.nrRows();
  }

  if (transb == 'N') {
    assert(b.nrRows() == k);
    assert(b.nrCols() == n);
  }
  else {
    assert(b.nrCols() == k);
    assert(b.nrRows() == n);
  }

  int lda = a.leadingDimension();
  int ldb = b.leadingDimension();
  int ldc = c.leadingDimension();

  blas::UseDevice<device_name>::gemm(&transa, &transb, m, n, k, alpha, a.ptr(), lda, b.ptr(), ldb,
                                     beta, c.ptr(), ldc, thread_id, stream_id);
}

// Performs the matrix-matrix multiplication c <- a * b
// Out: c
// Preconditions: a.nrRows() == c.nrRows(), b.nrCols() == c.nrCols() and a.nrCols() == b.nrRows()
template <typename Scalar, DeviceType device_name, template <typename, DeviceType> class MatrixA,
          template <typename, DeviceType> class MatrixB, template <typename, DeviceType> class MatrixC>
inline void gemm(const MatrixA<Scalar, device_name>& a, const MatrixB<Scalar, device_name>& b,
                 MatrixC<Scalar, device_name>& c, int thread_id = 0, int stream_id = 0) {
  gemm<Scalar, device_name>('N', 'N', 1., a, b, 0., c, thread_id, stream_id);
}

// Performs the matrix-matrix multiplication c <- alpha * a * b + beta * c,
// In/Out: c ('In' only if beta != 0)
// Preconditions: a.nrRows() == c.nrRows(), b.nrCols() == c.nrCols() and a.nrCols() == b.nrRows()
template <typename Scalar, DeviceType device_name, template <typename, DeviceType> class MatrixA,
          template <typename, DeviceType> class MatrixB, template <typename, DeviceType> class MatrixC>
inline void gemm(Scalar alpha, const MatrixA<Scalar, device_name>& a,
                 const MatrixB<Scalar, device_name>& b, Scalar beta,
                 MatrixC<Scalar, device_name>& c, int thread_id = 0, int stream_id = 0) {
  gemm<Scalar, device_name>('N', 'N', alpha, a, b, beta, c, thread_id, stream_id);
}

// Performs the matrix-matrix multiplication c <- op(a) * op(b),
// where op(X) = X if transX == 'N', op(X) = transposed(X) if transX == 'T', and
// op(X) == conjugate_transposed(X) if transX == 'C' (X = a, b).
// Out: c
// Preconditions: transa and transb should be one of the following: 'N', 'T' or 'C',
//                a.nrRows() == c.nrRows() if transa == 'N', a.nrCols() == c.nrRows() otherwise,
//                b.nrCols() == c.nrCols() if transb == 'N', b.nrRows() == c.nrCols() otherwise,
//                ka == kb, where ka = a.nrCols() if transa == 'N', ka = a.nrRows() otherwise and
//                          kb = b.nrRows() if transb == 'N', kb = b.nrCols() otherwise.
template <typename Scalar, DeviceType device_name>
inline void gemm(char transa, char transb, const Matrix<Scalar, device_name>& a,
                 const Matrix<Scalar, device_name>& b, Matrix<Scalar, device_name>& c,
                 int thread_id = 0, int stream_id = 0) {
  gemm<Scalar, device_name>(transa, transb, 1., a, b, 0., c, thread_id, stream_id);
}

// Performs the triangular solve b <- a^-1 * b,
// where a is a lower triangular matrix (uplo = 'L') or an upper triangular matrix (uplo = 'U'),
// with unit diagonal (diag = "U") or with general diagonal (diag = "N")
// In/Out: b
// Preconditions: a.nrRows() == a.nrCols() , a.nrCols() == b.nrRows()
template <typename Scalar, DeviceType device_name>
void trsm(char uplo, char diag, const Matrix<Scalar, device_name>& a,
          Matrix<Scalar, device_name>& b, int thread_id = 0, int stream_id = 0) {
  assert(uplo == 'U' or uplo == 'L');
  assert(diag == 'U' or diag == 'N');
  assert(a.nrRows() == a.nrCols());
  assert(b.nrRows() == a.nrCols());

  blas::UseDevice<device_name>::trsm("L", &uplo, "N", &diag, b.nrRows(), b.nrCols(), Scalar(1),
                                     a.ptr(), a.leadingDimension(), b.ptr(), b.leadingDimension(),
                                     thread_id, stream_id);
}

// Mixed real and complex matrix-matrix multiply.
// TODO: Not sure if this are needed.
//       The only file in which it may be used is basis_transformation_cd_to_ed.h.

// Performs the matrix-matrix multiplication c <- op(a) * op(b),
// where op(X) = X if transX == 'N', op(X) = transposed(X) if transX == 'T', and
// op(X) == conjugate_transposed(X) if transX == 'C' (X = a, b).
// Out: c
// Preconditions: transa should be one of the following: 'N', 'T',
//                transb should be one of 'N', 'T', 'C',
//                a.nrRows() == c.nrRows() if transa == 'N', a.nrCols() == c.nrRows() otherwise,
//                b.nrCols() == c.nrCols() if transb == 'N', b.nrRows() == c.nrCols() otherwise,
//                ka == kb, where ka = a.nrCols() if transa == 'N', ka = a.nrRows() otherwise and
//                          kb = b.nrRows() if transb == 'N', kb = b.nrCols() otherwise.
template <typename Scalar>
void gemm(char transa, char transb, Matrix<Scalar, CPU>& a, Matrix<std::complex<Scalar>, CPU>& b,
          Matrix<std::complex<Scalar>, CPU>& c) {
  Matrix<Scalar, CPU> b_part(b.size());
  Matrix<Scalar, CPU> c_re(c.size());
  Matrix<Scalar, CPU> c_im(c.size());

  Scalar sign = 1;
  if (transb == 'C') {
    sign = -1;
    transb = 'T';
  }

  for (int j = 0; j < b.nrCols(); ++j)
    for (int i = 0; i < b.nrRows(); ++i)
      b_part(i, j) = b(i, j).real();

  gemm(transa, transb, a, b_part, c_re);

  for (int j = 0; j < b.nrCols(); ++j)
    for (int i = 0; i < b.nrRows(); ++i)
      b_part(i, j) = b(i, j).imag();

  gemm(transa, transb, sign, a, b_part, Scalar(0), c_im);

  for (int j = 0; j < c.nrCols(); ++j)
    for (int i = 0; i < c.nrRows(); ++i)
      c(i, j) = std::complex<Scalar>(c_re(i, j), c_im(i, j));
}

// Performs the matrix-matrix multiplication c <- op(a) * op(b),
// where op(X) = X if transX == 'N', op(X) = transposed(X) if transX == 'T', and
// op(X) == conjugate_transposed(X) if transX == 'C' (X = a, b).
// Out: c
// Preconditions: transa should be one of the following: 'N', 'T', 'C',
//                transb should be one of 'N', 'T',
//                a.nrRows() == c.nrRows() if transa == 'N', a.nrCols() == c.nrRows() otherwise,
//                b.nrCols() == c.nrCols() if transb == 'N', b.nrRows() == c.nrCols() otherwise,
//                ka == kb, where ka = a.nrCols() if transa == 'N', ka = a.nrRows() otherwise and
//                          kb = b.nrRows() if transb == 'N', kb = b.nrCols() otherwise.
template <typename Scalar>
static void gemm(char transa, char transb, Matrix<std::complex<Scalar>, CPU>& a,
                 Matrix<Scalar, CPU>& b, Matrix<std::complex<Scalar>, CPU>& c) {
  Matrix<Scalar, CPU> a_part(a.size());
  Matrix<Scalar, CPU> c_re(c.size());
  Matrix<Scalar, CPU> c_im(c.size());

  Scalar sign = 1;
  if (transa == 'C') {
    sign = -1;
    transa = 'T';
  }

  for (int j = 0; j < a.nrCols(); ++j)
    for (int i = 0; i < a.nrRows(); ++i)
      a_part(i, j) = a(i, j).real();

  gemm(transa, transb, a_part, b, c_re);

  for (int j = 0; j < a.nrCols(); ++j)
    for (int i = 0; i < a.nrRows(); ++i)
      a_part(i, j) = a(i, j).imag();

  gemm(transa, transb, sign, a_part, b, Scalar(0), c_im);

  for (int j = 0; j < c.nrCols(); ++j)
    for (int i = 0; i < c.nrRows(); ++i)
      c(i, j) = std::complex<Scalar>(c_re(i, j), c_im(i, j));
}

// Performs the matrix-matrix multiplication c = op(a) * op(b), where each matrix is split in real
// and imaginary part. This is implemented with 3 real matrix-matrix multiplications.
// Out: c
// Preconditions: transa and transb should be one of the following: 'N', 'T', 'C'.
//                a[0].size == a[1].size()
//                b[0].size == b[1].size()
//                c[0].size == c[1].size()
//                a[0].nrRows() == c[0].nrRows() if transa == 'N', a[0].nrCols() == c[0].nrRows()
//                  otherwise,
//                b[0].nrCols() == c[0].nrCols() if transb == 'N', b[0].nrRows() == c[0].nrCols()
//                  otherwise,
//                ka == kb, where ka = a[0].nrCols() if transa == 'N', ka = a[0].nrRows() otherwise
//                and kb = b[0].nrRows() if transb == 'N', kb = b[0].nrCols() otherwise.
template <typename Scalar>
void multiply(char transa, char transb, const std::array<Matrix<Scalar, CPU>, 2>& a,
              const std::array<Matrix<Scalar, CPU>, 2>& b, std::array<Matrix<Scalar, CPU>, 2>& c,
              std::array<Matrix<Scalar, CPU>, 5>& work) {
  assert(a[0].size() == a[1].size());
  assert(b[0].size() == b[1].size());
  assert(c[0].size() == c[1].size());

  work[0].resizeNoCopy(c[0].size());
  work[1].resizeNoCopy(c[0].size());
  work[2].resizeNoCopy(c[0].size());
  auto& a_sum = work[3];
  auto& b_sum = work[4];
  a_sum.resizeNoCopy(a[0].size());
  b_sum.resizeNoCopy(b[0].size());

  const Scalar signa = transa == 'C' ? transa = 'T', -1 : 1;
  const Scalar signb = transb == 'C' ? transb = 'T', -1 : 1;

  for (int j = 0; j < a[0].nrCols(); ++j)
    for (int i = 0; i < a[0].nrRows(); ++i)
      a_sum(i, j) = a[0](i, j) + signa * a[1](i, j);
  for (int j = 0; j < b[0].nrCols(); ++j)
    for (int i = 0; i < b[0].nrRows(); ++i)
      b_sum(i, j) = b[0](i, j) + signb * b[1](i, j);

  gemm(transa, transb, a[0], b[0], work[0]);
  gemm(transa, transb, signa * signb, a[1], b[1], Scalar(0), work[1]);
  gemm(transa, transb, a_sum, b_sum, work[2]);

  for (int j = 0; j < c[0].nrCols(); ++j)
    for (int i = 0; i < c[0].nrRows(); ++i) {
      c[0](i, j) = work[0](i, j) - work[1](i, j);
      c[1](i, j) = work[2](i, j) - work[0](i, j) - work[1](i, j);
    }
}

template <typename Scalar>
void multiply(const std::array<Matrix<Scalar, CPU>, 2>& a,
              const std::array<Matrix<Scalar, CPU>, 2>& b, std::array<Matrix<Scalar, CPU>, 2>& c,
              std::array<Matrix<Scalar, CPU>, 5>& work) {
  multiply('N', 'N', a, b, c, work);
}

// Performs the matrix-matrix multiplication c = op(a) * op(b), where a and c are split in real
// and imaginary part, while b is real.
// Out: c
// Preconditions: transa and transb should be one of the following: 'N', 'T',
//                a[0].size == a[1].size()
//                c[0].size == c[1].size()
//                a[0].nrRows() == c[0].nrRows() if transa == 'N', a.nrCols() == c.nrRows()
//                  otherwise,
//                b.nrCols() == c[0].nrCols() if transb == 'N', b.nrRows() == c.nrCols()
//                  otherwise,
//                ka == kb, where ka = a[0].nrCols() if transa == 'N', ka = a[0].nrRows() otherwise
//                and kb = b.nrRows() if transb == 'N', kb = b.nrCols() otherwise.
template <typename Scalar, DeviceType device_name>
void multiply(char transa, char transb, const std::array<Matrix<Scalar, device_name>, 2>& a,
              const Matrix<Scalar, device_name>& b, std::array<Matrix<Scalar, device_name>, 2>& c) {
  assert(transa == 'N' || transa == 'T' || transa == 'C');
  assert(transb == 'N' || transb == 'T');
  assert(a[0].size() == a[1].size());
  assert(c[0].size() == c[1].size());

  gemm(transa, transb, a[0], b, c[0]);
  const Scalar sign = transa == 'C' ? transa = 'T', -1 : 1;
  gemm(transa, transb, sign, a[1], b, Scalar(0), c[1]);
}

template <typename Scalar, DeviceType device_name>
void multiply(const std::array<Matrix<Scalar, device_name>, 2>& a,
              const Matrix<Scalar, device_name>& b, std::array<Matrix<Scalar, device_name>, 2>& c) {
  multiply('N', 'N', a, b, c);
}

// Performs the matrix-matrix multiplication c = op(a) * op(b), where b and c are split in real
// and imaginary part, while a is real.
// Out: c
// Preconditions: transa and transb should be one of the following: 'N', 'T',
//                b[0].size == b[1].size()
//                c[0].size == c[1].size()
//                a.nrRows() == c[0].nrRows() if transa == 'N', a.nrCols() == c.nrRows()
//                  otherwise,
//                b.[0]nrCols() == c[0].nrCols() if transb == 'N', b.[0]nrRows() == c.nrCols()
//                  otherwise,
//                ka == kb, where ka = a.nrCols() if transa == 'N', ka = a.nrRows() otherwise
//                and kb = b.[0]nrRows() if transb == 'N', kb = b.[0]nrCols() otherwise.
template <typename Scalar, DeviceType device_name>
void multiply(char transa, char transb, const Matrix<Scalar, device_name>& a,
              const std::array<Matrix<Scalar, device_name>, 2>& b,
              std::array<Matrix<Scalar, device_name>, 2>& c) {
  assert(transa == 'N' || transa == 'T');
  assert(transb == 'N' || transb == 'T' || transb == 'C');
  assert(b[0].size() == b[1].size());
  assert(c[0].size() == c[1].size());

  gemm(transa, transb, a, b[0], c[0]);
  const Scalar sign = transb == 'C' ? transb = 'T', -1 : 1;
  gemm(transa, transb, sign, a, b[1], Scalar(0), c[1]);
}

template <typename Scalar, DeviceType device_name>
void multiply(const Matrix<Scalar, device_name>& a,
              const std::array<Matrix<Scalar, device_name>, 2>& b,
              std::array<Matrix<Scalar, device_name>, 2>& c) {
  multiply('N', 'N', a, b, c);
}

// Performs the matrix-matrix multiplication b <- D * a,
// where d is a vector containing the diagonal elements of the matrix D.
// Out: b
// Preconditions: a.size() == b.size(), d.size() == a.nrRows().
template <typename ScalarIn, typename ScalarOut, DeviceType device_name>
inline void multiplyDiagonalLeft(const Vector<ScalarIn, device_name>& d,
                                 const Matrix<ScalarIn, device_name>& a,
                                 Matrix<ScalarOut, device_name>& b, int thread_id = 0,
                                 int stream_id = 0) {
  lapack::UseDevice<device_name>::multiplyDiagonalLeft(a.nrRows(), a.nrCols(), d.ptr(), 1, a.ptr(),
                                                       a.leadingDimension(), b.ptr(),
                                                       b.leadingDimension(), thread_id, stream_id);
}
template <typename ScalarIn, typename ScalarOut>
inline void multiplyDiagonalLeft(const Vector<ScalarIn, CPU>& d, const Matrix<ScalarIn, GPU>& a,
                                 Matrix<ScalarOut, GPU>& b, int thread_id = 0, int stream_id = 0) {
  Vector<ScalarIn, GPU> d_gpu(d);
  multiplyDiagonalLeft(d_gpu, a, b, thread_id, stream_id);
}

// Performs the matrix-matrix multiplication b <- a * D,
// where d is a vector containing the diagonal elements of the matrix D.
// Out: b
// Preconditions: a.size() == b.size(), d.size() == a.nrCols().
template <typename Scalar, DeviceType device_name>
inline void multiplyDiagonalRight(const Matrix<Scalar, device_name>& a,
                                  const Vector<Scalar, device_name>& d,
                                  Matrix<Scalar, device_name>& b, int thread_id = 0,
                                  int stream_id = 0) {
  lapack::UseDevice<device_name>::multiplyDiagonalRight(a.nrRows(), a.nrCols(), a.ptr(),
                                                        a.leadingDimension(), d.ptr(), 1, b.ptr(),
                                                        b.leadingDimension(), thread_id, stream_id);
}
template <typename Scalar>
inline void multiplyDiagonalRight(const Matrix<Scalar, GPU>& a, const Vector<Scalar, CPU>& d,
                                  Matrix<Scalar, GPU>& b, int thread_id = 0, int stream_id = 0) {
  Vector<Scalar, GPU> d_gpu(d);
  multiplyDiagonalRight(a, d_gpu, b, thread_id, stream_id);
}

// Computes the eigenvalues, the left eigenvectors (if jobvl == 'V')
// and the right eigenvectors (if jobvr == 'V') of the real matrix a.
// The real parts of the eigenvalues are stored in lambda_re, while the imaginary parts in
// lambda_im.
// If computed the left eigenvectors are stored in vl and the right eigenvectors in vr.
// See sgeev, dgeev Lapack documentation for information about how the
// eigenvectors are stored.
// Out: lambda_re, lambda_im, vl, vr.
// Precondition: jobvl == 'N' or jobvl == 'V',
//               jobvr == 'N' or jobvr == 'V',
//               a is a square matrix.
// Postcondition: lambda_re, lambda_i, are resized, vl if jobvl == 'V', vr if jobvr == 'V' are
// resized.
template <typename Scalar>
void eigensolver(char jobvl, char jobvr, const Matrix<Scalar, CPU>& a, Vector<Scalar, CPU>& lambda_re,
                 Vector<Scalar, CPU>& lambda_im, Matrix<Scalar, CPU>& vl, Matrix<Scalar, CPU>& vr) {
  assert(a.is_square());

  Matrix<Scalar, CPU> a_copy(a);
  lambda_re.resizeNoCopy(a_copy.nrRows());
  lambda_im.resizeNoCopy(a_copy.nrRows());
  int ldvl = 1;
  int ldvr = 1;
  if (jobvl == 'V' || jobvl == 'v') {
    vl.resizeNoCopy(a_copy.size());
    ldvl = vl.leadingDimension();
  }
  if (jobvr == 'V' || jobvr == 'v') {
    vr.resizeNoCopy(a_copy.size());
    ldvr = vr.leadingDimension();
  }

  // Get optimal worksize.
  int lwork = util::getEigensolverWorkSize(jobvl, jobvr, a_copy);
  dca::linalg::Vector<Scalar, CPU> work(lwork);

  lapack::geev(&jobvl, &jobvr, a_copy.nrRows(), a_copy.ptr(), a_copy.leadingDimension(),
               lambda_re.ptr(), lambda_im.ptr(), vl.ptr(), ldvl, vr.ptr(), ldvr, work.ptr(),
               work.size());
}

// Computes the eigenvalues, the left eigenvectors (if jobvl == 'V')
// and the right eigenvectors (if jobvr == 'V') of the complex matrix a.
// The eigenvalues are stored in lambda.
// If computed the left eigenvectors are stored in vl and the right eigenvectors in vr.
// Out: lambda, vl, vr.
// Precondition: jobvl == 'N' or jobvl == 'V',
//               jobvr == 'N' or jobvr == 'V',
//               a is a square matrix.
// Postcondition: lambda, is resized, vl if jobvl == 'V', vr if jobvr == 'V' are resized.
template <typename Scalar>
void eigensolver(char jobvl, char jobvr, const Matrix<std::complex<Scalar>, CPU>& a,
                 Vector<std::complex<Scalar>, CPU>& lambda, Matrix<std::complex<Scalar>, CPU>& vl,
                 Matrix<std::complex<Scalar>, CPU>& vr) {
  assert(a.is_square());

  Matrix<std::complex<Scalar>, CPU> a_copy(a);
  lambda.resizeNoCopy(a_copy.nrRows());
  int ldvl = 1;
  int ldvr = 1;
  if (jobvl == 'V' || jobvl == 'v') {
    vl.resizeNoCopy(a_copy.size());
    ldvl = vl.leadingDimension();
  }
  if (jobvr == 'V' || jobvr == 'v') {
    vr.resizeNoCopy(a_copy.size());
    ldvr = vr.leadingDimension();
  }

  // Get optimal worksize.
  int lwork = util::getEigensolverWorkSize(jobvl, jobvr, a_copy);
  dca::linalg::Vector<std::complex<Scalar>, CPU> work(lwork);
  dca::linalg::Vector<Scalar, CPU> rwork(2 * a_copy.nrRows());

  lapack::geev(&jobvl, &jobvr, a_copy.nrRows(), a_copy.ptr(), a_copy.leadingDimension(),
               lambda.ptr(), vl.ptr(), ldvl, vr.ptr(), ldvr, work.ptr(), work.size(), rwork.ptr());
}

// Computes the eigenvalues, and the eigenvectors (if jobv == 'V') of the real symmetric matrix a.
// if uplo == 'U' the upper triangular part of a is referenced, whereas
// if uplo == 'L' the lower triangular part of a is referenced.
// The eigenvalues are stored in lambda.
// If computed the eigenvectors are stored in v.
// Out: lambda, v
// Precondition: jobv == 'N' or jobv == 'V',
//               uplo == 'U' or uplo == 'L',
//               a is a square matrix.
// Postcondition: lambda, and v are resized.
template <typename Scalar>
void eigensolverSymmetric(char jobv, char uplo, const Matrix<Scalar, CPU>& a,
                          Vector<Scalar, CPU>& lambda, Matrix<Scalar, CPU>& v) {
  assert(a.is_square());

  lambda.resizeNoCopy(a.nrRows());
  v = a;

  // Get optimal worksize.
  auto lwork = util::getEigensolverSymmetricWorkSize(jobv, uplo, v);
  dca::linalg::Vector<Scalar, CPU> work(std::get<0>(lwork));
  dca::linalg::Vector<int, CPU> iwork(std::get<1>(lwork));

  lapack::syevd(&jobv, &uplo, v.nrRows(), v.ptr(), v.leadingDimension(), lambda.ptr(), work.ptr(),
                work.size(), iwork.ptr(), iwork.size());
}
// For real types Hermitian and symmetric is the same.
template <typename Scalar>
inline void eigensolverHermitian(char jobv, char uplo, const Matrix<Scalar, CPU>& a,
                                 Vector<Scalar, CPU>& lambda, Matrix<Scalar, CPU>& v) {
  eigensolverSymmetric(jobv, uplo, a, lambda, v);
}

// Computes the eigenvalues, and the eigenvectors (if jobv == 'V')
// of the complex Hermitian matrix a.
// if uplo == 'U' the upper triangular part of a is referenced, whereas
// if uplo == 'L' the lower triangular part of a is referenced.
// The eigenvalues are stored in lambda.
// If computed the eigenvectors are stored in v.
// Out: lambda, v
// Precondition: jobv == 'N' or jobv == 'V',
//               uplo == 'U' or uplo == 'L',
//               a is a square matrix.
// Postcondition: lambda, and v are resized.
template <typename Scalar>
void eigensolverHermitian(char jobv, char uplo, const Matrix<std::complex<Scalar>, CPU>& a,
                          Vector<Scalar, CPU>& lambda, Matrix<std::complex<Scalar>, CPU>& v) {
  assert(a.is_square());

  lambda.resizeNoCopy(a.nrRows());
  v = a;

  // Get optimal worksize.
  auto lwork = util::getEigensolverHermitianWorkSize(jobv, uplo, v);
  dca::linalg::Vector<std::complex<Scalar>, CPU> work(std::get<0>(lwork));
  dca::linalg::Vector<Scalar, CPU> rwork(std::get<1>(lwork));
  dca::linalg::Vector<int, CPU> iwork(std::get<2>(lwork));

  lapack::heevd(&jobv, &uplo, v.nrRows(), v.ptr(), v.leadingDimension(), lambda.ptr(), work.ptr(),
                work.size(), rwork.ptr(), rwork.size(), iwork.ptr(), iwork.size());
}

template <typename Scalar>
void eigensolverGreensFunctionMatrix(char jobv, char uplo, const Matrix<std::complex<Scalar>, CPU>& a,
                                     Vector<Scalar, CPU>& lambda,
                                     Matrix<std::complex<Scalar>, CPU>& v) {
  assert(a.is_square());
  int n = a.nrRows();
  assert(n % 2 == 0);

  if (n == 2) {
    lambda.resize(2);
    v.resize(2);

    assert(std::abs(std::imag(a(0, 0))) < 1.e-6);
    assert(std::abs(std::imag(a(1, 1))) < 1.e-6);

    lambda[0] = std::real(a(0, 0));
    lambda[1] = std::real(a(1, 1));
    v(0, 0) = 1.;
    v(1, 0) = 0.;
    v(0, 1) = 0.;
    v(1, 1) = 1.;
    return;
  }
  eigensolverHermitian(jobv, uplo, a, lambda, v);
}

// Computes the pseudo inverse of the matrix.
// Out: a_inv
// Postconditions: a_inv is resized to the needed dimension.
template <typename Scalar>
void pseudoInverse(const Matrix<Scalar, CPU>& a, Matrix<Scalar, CPU>& a_inv, double eps = 1.e-6) {
  int m = a.nrRows();
  int n = a.nrCols();
  a_inv.resizeNoCopy(std::make_pair(n, m));

  using RealType = decltype(std::real(*a.ptr()));

  if (m <= n) {
    // a_inv = a'*inv(a*a')
    // inv(a*a') = v*inv(lambda)*v', [lambda, v] = eig(a*a')

    Matrix<Scalar, CPU> a_at("A_At", m);
    dca::linalg::matrixop::gemm('N', 'C', a, a, a_at);

    dca::linalg::Vector<RealType, CPU> lambda("Lambda", m);
    Matrix<Scalar, CPU> v("V", m);

    eigensolverHermitian('V', 'U', a_at, lambda, v);
    Matrix<Scalar, CPU> vt(v);

    for (int j = 0; j < m; j++) {
      Scalar lambda_inv = 0;

      if (lambda[j] > eps * lambda[m - 1])
        lambda_inv = 1. / lambda[j];

      scaleCol(v, j, lambda_inv);
    }

    gemm('N', 'C', v, vt, a_at);
    gemm('C', 'N', a, a_at, a_inv);
  }
  else {
    // a_inv = inv(a'*a)*a'
    // inv(a'*a) = v*inv(lambda)*v', [lambda, v] = eig(a'*a)

    Matrix<Scalar, CPU> at_a("at_a", n);
    dca::linalg::matrixop::gemm('C', 'N', a, a, at_a);

    dca::linalg::Vector<RealType, CPU> lambda("Lambda", n);
    Matrix<Scalar, CPU> v("V", n);

    eigensolverHermitian('V', 'U', at_a, lambda, v);
    Matrix<Scalar, CPU> vt(v);

    for (int j = 0; j < n; j++) {
      Scalar lambda_inv = 0;

      if (lambda[j] > eps * lambda[n - 1])
        lambda_inv = 1. / lambda[j];

      scaleCol(v, j, lambda_inv);
    }

    gemm('N', 'C', v, vt, at_a);
    gemm('N', 'C', at_a, a, a_inv);
  }
}

// Computes (in place) the determinant of the matrix.
// Returns: determinant.
// Postcondition: M is its LU decomposition.
template <template <typename, DeviceType> class MatrixType, typename Scalar>
Scalar determinantIP(MatrixType<Scalar, CPU>& M) {
  assert(M.nrCols() == M.nrRows());
  const int n = M.nrCols();
  std::vector<int> ipiv(n);

  try {
    lapack::getrf(n, n, M.ptr(), M.leadingDimension(), ipiv.data());
  }
  catch (lapack::util::LapackException& err) {
    if (err.info() > 0)
      return 0;
    else
      throw(std::logic_error("LU decomposition failed."));
  }

  double det = 1.;
  for (int i = 0; i < n; i++) {
    det *= M(i, i);
    if (ipiv[i] != i + 1)
      det *= -1;
  }
  return det;
}

// Copy and computes the determinant of the matrix.
// Returns: determinant.
template <typename Scalar, DeviceType device>
double determinant(const Matrix<Scalar, device>& M) {
  Matrix<Scalar, CPU> M_copy(M);
  return determinantIP(M_copy);
}

// Returns: logarithm of the absolute value of the determinant and the sign of the determinant,
//          or zero if the determinant is zero.
// Postcondition: M is its LU decomposition.
template <template <typename, DeviceType> class MatrixType, typename Scalar>
std::pair<Scalar, int> logDeterminantIP(MatrixType<Scalar, CPU>& M, std::vector<int>& ipiv) {
  assert(M.is_square());
  static_assert(std::is_same_v<Scalar, float> || std::is_same_v<Scalar, double>,
                " This function is defined only for Real numbers");

  const int n = M.nrCols();
  ipiv.resize(n);

  try {
    lapack::getrf(n, n, M.ptr(), M.leadingDimension(), ipiv.data());
  }
  catch (lapack::util::LapackException& err) {
    if (err.info() > 0)
      return {0., 0};
    else
      throw(std::logic_error("LU decomposition failed."));
  }

  Scalar log_det = 0.;
  int sign = 1;

  for (int i = 0; i < n; i++) {
    log_det += std::log(std::abs(M(i, i)));
    if (M(i, i) < 0)
      sign *= -1;

    if (ipiv[i] != i + 1)
      sign *= -1;
  }

  return {log_det, sign};
}

template <template <typename, DeviceType> class MatrixType, typename Scalar, DeviceType device>
auto logDeterminant(const MatrixType<Scalar, device>& m) {
  Matrix<Scalar, CPU> m_copy(m);
  std::vector<int> ipiv;
  return logDeterminantIP(m_copy, ipiv);
}

template <typename Scalar, template <typename, DeviceType> class MatrixType>
std::pair<Scalar, int> inverseAndLogDeterminant(MatrixType<Scalar, CPU>& mat) {
  std::vector<int> ipiv;
  const auto [log_det, sign] = logDeterminantIP(mat, ipiv);

  if (!sign)
    throw(std::logic_error("Singular matrix"));

  const int lwork = util::getInverseWorkSize(mat);
  std::vector<Scalar> work(lwork);

  lapack::UseDevice<CPU>::getri(mat.nrRows(), mat.ptr(), mat.leadingDimension(), ipiv.data(),
                                work.data(), lwork);

  return {-log_det, sign};
}

template <typename Scalar>
bool areNear(const Matrix<Scalar, CPU>& A, const Matrix<Scalar, CPU>& B, const double err = 1e-16) {
  if (A.size() != B.size())
    return false;

  for (int j = 0; j < A.size().second; j++)
    for (int i = 0; i < A.size().first; i++)
      if (std::abs(A(i, j) - B(i, j)) > err)
        return false;

  return true;
}

}  // namespace matrixop
}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_MATRIXOP_HPP
