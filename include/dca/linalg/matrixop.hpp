// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the matrix interface for the following matrix operations:
// - copyCol, copyRow, copyCols, copyRows
// - difference
// - insertCol, insertRow (for CPU matrices only)
// - inverse
// - removeCol, removeCols, removeRow, removeRows, removeRowAndCol, removeRowAndCols
// - scaleCol, scaleRow, scaleRows
// - swapCol, swapRow, swapRowAndCol
// - swapCols, swapRows (for GPU matrices only)
// - gemm
// - multiply
// - trsm
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
template <typename ScalarType>
inline void copyMatrixToArray(const Matrix<ScalarType, CPU>& mat, ScalarType* a, int lda) {
  assert(lda >= mat.nrRows());
  lapack::lacpy("A", mat.nrRows(), mat.nrCols(), mat.ptr(), mat.leadingDimension(), a, lda);
}

// Copies the m by n matrix stored in a to the matrix mat.
// Preconditions: lda >= m.
template <typename ScalarType>
inline void copyArrayToMatrix(int m, int n, const ScalarType* a, int lda,
                              Matrix<ScalarType, CPU>& mat) {
  assert(lda >= m);
  mat.resizeNoCopy(std::make_pair(m, n));
  lapack::lacpy("A", mat.nrRows(), mat.nrCols(), a, lda, mat.ptr(), mat.leadingDimension());
}

// Copies the jx-th column of mat_x into the jy-th column of mat_y.
// In/Out: mat_y
// Preconditions: mat_x.nrRows() == mat_y.nrRows(),
//                0 <= jx < mat_x.nrCols(), 0 <= jy < mat_y.nrCols().
template <typename ScalarType, DeviceType device_name>
inline void copyCol(const Matrix<ScalarType, device_name>& mat_x, int jx,
                    Matrix<ScalarType, device_name>& mat_y, int jy, int thread_id = 0,
                    int stream_id = 0) {
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
template <typename ScalarType>
inline void copyCols(const Matrix<ScalarType, CPU>& mat_x, const Vector<int, CPU>& j_x,
                     Matrix<ScalarType, CPU>& mat_y, const Vector<int, CPU>& j_y,
                     int /*thread_id*/ = 0, int /*stream_id*/ = 0) {
  assert(j_x.size() <= j_y.size());

  for (int ind_j = 0; ind_j < j_x.size(); ++ind_j)
    copyCol(mat_x, j_x[ind_j], mat_y, j_y[ind_j]);
}
#ifdef DCA_HAVE_CUDA
template <typename ScalarType>
inline void copyCols(const Matrix<ScalarType, GPU>& mat_x, const Vector<int, GPU>& j_x,
                     Matrix<ScalarType, GPU>& mat_y, const Vector<int, GPU>& j_y, int thread_id = 0,
                     int stream_id = 0) {
  assert(j_x.size() <= j_y.size());
  assert(mat_x.nrRows() == mat_y.nrRows());

  blas::copyCols(mat_x.nrRows(), j_x.size(), j_x.ptr(), mat_x.ptr(), mat_x.leadingDimension(),
                 j_y.ptr(), mat_y.ptr(), mat_y.leadingDimension(), thread_id, stream_id);
}
#endif  // DCA_HAVE_CUDA

// Copies the ix-th row of mat_x into the iy-th row of mat_y.
// In/Out: mat_y
// Preconditions: mat_x.nrCols() == mat_y.nrCols(),
//                0 <= ix < mat_x.nrRows(), 0 <= iy < mat_y.nrRows().
template <typename ScalarType, DeviceType device_name>
inline void copyRow(const Matrix<ScalarType, device_name>& mat_x, int ix,
                    Matrix<ScalarType, device_name>& mat_y, int iy, int thread_id = 0,
                    int stream_id = 0) {
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
template <typename ScalarType>
inline void copyRows(const Matrix<ScalarType, CPU>& mat_x, const Vector<int, CPU>& i_x,
                     Matrix<ScalarType, CPU>& mat_y, const Vector<int, CPU>& i_y,
                     int /*thread_id*/ = 0, int /*stream_id*/ = 0) {
  assert(i_x.size() <= i_y.size());
  assert(mat_x.nrCols() == mat_y.nrCols());

  for (int j = 0; j < mat_x.nrCols(); ++j)
    for (int ind_i = 0; ind_i < i_x.size(); ++ind_i)
      mat_y(i_y[ind_i], j) = mat_x(i_x[ind_i], j);
}
#ifdef DCA_HAVE_CUDA
template <typename ScalarType>
inline void copyRows(const Matrix<ScalarType, GPU>& mat_x, const Vector<int, GPU>& i_x,
                     Matrix<ScalarType, GPU>& mat_y, const Vector<int, GPU>& i_y, int thread_id = 0,
                     int stream_id = 0) {
  assert(i_x.size() <= i_y.size());
  assert(mat_x.nrCols() == mat_y.nrCols());

  blas::copyRows(mat_x.nrCols(), i_x.size(), i_x.ptr(), mat_x.ptr(), mat_x.leadingDimension(),
                 i_y.ptr(), mat_y.ptr(), mat_y.leadingDimension(), thread_id, stream_id);
}
#endif  // DCA_HAVE_CUDA

// Returns the difference of two matrices in terms of max_i,j(|a(i, j) - b(i, j)|).
// If the difference is larger than the threshold a std::logic_error exception is thrown,
// and if NDEBUG is not defined each difference which exceeds the threshold is printed.
// Preconditions: a.size() == b.size().
template <typename ScalarType>
auto difference(const Matrix<ScalarType, CPU>& a, const Matrix<ScalarType, CPU>& b,
                double diff_threshold = 1e-3) {
  assert(a.size() == b.size());

  auto max_diff = std::abs(ScalarType(0));

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
template <typename ScalarType, DeviceType device_name>
auto difference(const Matrix<ScalarType, device_name>& a, const Matrix<ScalarType, CPU>& b,
                double diff_threshold = 1e-3) {
  Matrix<ScalarType, CPU> cp_a(a);
  return difference(cp_a, b, diff_threshold);
}
template <typename ScalarType, DeviceType device_name_a, DeviceType device_name_b>
auto difference(const Matrix<ScalarType, device_name_a>& a,
                const Matrix<ScalarType, device_name_b>& b, double diff_threshold = 1e-3) {
  Matrix<ScalarType, CPU> cp_b(b);
  return difference(a, cp_b, diff_threshold);
}

// Insert a column at position j. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= j < mat.nrCols() + 1.
// Postconditions: The elements of the inserted column are set to 0.
template <typename ScalarType>
void insertCol(Matrix<ScalarType, CPU>& mat, int j) {
  assert(j >= 0 && j < mat.nrCols() + 1);

  mat.resize(std::make_pair(mat.nrRows(), mat.nrCols() + 1));

  if (mat.nrRows() > 0 && j < mat.nrCols() - 1)
    memmove(mat.ptr(0, j + 1), mat.ptr(0, j),
            sizeof(ScalarType) * (mat.nrCols() - 1 - j) * mat.leadingDimension());

  for (int i = 0; i < mat.nrRows(); ++i)
    mat(i, j) = 0;
}

// Insert a row at position i. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= i < mat.nrRows() + 1.
// Postconditions: The elements of the inserted row are set to 0.
template <typename ScalarType>
void insertRow(Matrix<ScalarType, CPU>& mat, int i) {
  assert(i >= 0 && i < mat.nrRows() + 1);

  mat.resize(std::make_pair(mat.nrRows() + 1, mat.nrCols()));

  if (i < mat.nrRows() - 1)
    for (int j = 0; j < mat.nrCols(); ++j)
      memmove(mat.ptr(i + 1, j), mat.ptr(i, j), sizeof(ScalarType) * (mat.nrRows() - 1 - i));

  for (int j = 0; j < mat.nrCols(); ++j)
    mat(i, j) = 0;
}

// Computes the inverse of the matrix using the LU factorization.
// In/Out: mat
// Out: ipiv, work
// Preconditions: mat is a square matrix.
// Postconditions: ipiv and work are resized to the needed dimension.
template <typename ScalarType, DeviceType device_name, template <typename, DeviceType> class MatrixType>
void inverse(MatrixType<ScalarType, device_name>& mat, Vector<int, CPU>& ipiv,
             Vector<ScalarType, device_name>& work) {
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
template <typename ScalarType, DeviceType device_name, template <typename, DeviceType> class MatrixType>
void inverse(MatrixType<ScalarType, device_name>& mat) {
  Vector<int, CPU> ipiv;
  Vector<ScalarType, device_name> work;
  inverse(mat, ipiv, work);
}

// Remove the j-th column. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= j < mat.nrCols().
template <typename ScalarType>
void removeCol(Matrix<ScalarType, CPU>& mat, int j) {
  assert(j >= 0 && j < mat.nrCols());

  if (mat.nrRows() > 0 && j < mat.nrCols() - 1)
    memmove(mat.ptr(0, j), mat.ptr(0, j + 1),
            sizeof(ScalarType) * (mat.nrCols() - j - 1) * mat.leadingDimension());

  mat.resize(std::make_pair(mat.nrRows(), mat.nrCols() - 1));
}

#ifdef DCA_HAVE_CUDA
template <typename ScalarType>
void removeCol(Matrix<ScalarType, GPU>& mat, int j) {
  assert(j >= 0 && j < mat.nrCols());

  if (mat.nrRows() > 0 && j < mat.nrCols() - 1)
    blas::moveLeft(mat.nrRows(), mat.nrCols() - j, mat.ptr(0, j), mat.leadingDimension());

  mat.resize(std::make_pair(mat.nrRows(), mat.nrCols() - 1));
}
#endif  // DCA_HAVE_CUDA

// Remove columns in range [first, last]. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= first, last < mat.nrCols().
template <typename ScalarType>
void removeCols(Matrix<ScalarType, CPU>& mat, int first, int last) {
  const int n_removed = last - first + 1;
  const int n = mat.nrRows();
  const int m = mat.nrCols();
  assert(last < m and last >= first and first >= 0);

  if (n > 0 and last < m - 1)
    std::memmove(mat.ptr(0, first), mat.ptr(0, last + 1),
                 mat.leadingDimension() * (m - last - 1) * sizeof(ScalarType));

  mat.resize(std::make_pair(n, m - n_removed));
}

// Remove the i-th row. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= i < mat.nrRows().
template <typename ScalarType>
void removeRow(Matrix<ScalarType, CPU>& mat, int i) {
  assert(i >= 0 && i < mat.nrRows());

  if (i < mat.nrRows() - 1)
    for (int j = 0; j < mat.nrCols(); ++j)
      memmove(mat.ptr(i, j), mat.ptr(i + 1, j), sizeof(ScalarType) * (mat.nrRows() - i - 1));

  mat.resize(std::make_pair(mat.nrRows() - 1, mat.nrCols()));
}

#ifdef DCA_HAVE_CUDA
template <typename ScalarType>
void removeRow(Matrix<ScalarType, GPU>& mat, int i) {
  assert(i >= 0 && i < mat.nrRows());

  if (mat.nrCols() > 0 && i < mat.nrRows() - 1)
    blas::moveUp(mat.nrRows() - i, mat.nrCols(), mat.ptr(i, 0), mat.leadingDimension());

  mat.resize(std::make_pair(mat.nrRows() - 1, mat.nrCols()));
}
#endif  // DCA_HAVE_CUDA

// Remove rows in range [first, last]. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= first, last < mat.nrRows().
template <typename ScalarType>
void removeRows(Matrix<ScalarType, CPU>& mat, int first, int last) {
  const int n_removed = last - first + 1;
  const int n = mat.nrRows();
  const int m = mat.nrCols();
  assert(last < n and last >= first and first >= 0);

  if (last < n - 1)
    for (int j = 0; j < m; ++j)
      std::memmove(mat.ptr(first, j), mat.ptr(last + 1, j), (n - last - 1) * sizeof(ScalarType));

  mat.resize(std::make_pair(n - n_removed, m));
}

// Remove the i-th row and the j-th column. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= i < mat.nrRows(), 0 <= j < mat.nrCols().
template <typename ScalarType, DeviceType device_name>
inline void removeRowAndCol(Matrix<ScalarType, device_name>& mat, int i, int j) {
  removeRow(mat, i);
  removeCol(mat, j);
}

// Remove the i-th row and the i-th column. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= i < mat.nrRows(), i < mat.nrCols().
template <typename ScalarType, DeviceType device_name>
inline void removeRowAndCol(Matrix<ScalarType, device_name>& mat, int i) {
  removeRowAndCol(mat, i, i);
}

// Remove rows and columns in range [first, last]. The data is moved accordingly.
// In/Out: mat
// Preconditions: 0 <= first, last < min(mat.nrRows(), mat.nrCols()).
template <typename ScalarType>
void removeRowsAndCols(Matrix<ScalarType, CPU>& mat, int first, int last) {
  removeCols(mat, first, last);
  removeRows(mat, first, last);
}

// Scales the j-th column of mat by val.
// In/Out: mat
// Preconditions: 0 <= j < mat.nrCols().
template <typename ScalarType, DeviceType device_name>
inline void scaleCol(Matrix<ScalarType, device_name>& mat, int j, ScalarType val, int thread_id = 0,
                     int stream_id = 0) {
  assert(j >= 0 && j < mat.nrCols());
  blas::UseDevice<device_name>::scal(mat.nrRows(), val, mat.ptr(0, j), 1, thread_id, stream_id);
}

// Scales the i-th row of mat by val.
// In/Out: mat
// Preconditions: 0 <= i < mat.nrRow().
template <typename ScalarType, DeviceType device_name>
inline void scaleRow(Matrix<ScalarType, device_name>& mat, int i, ScalarType val, int thread_id = 0,
                     int stream_id = 0) {
  assert(i >= 0 && i < mat.nrRows());
  blas::UseDevice<device_name>::scal(mat.nrCols(), val, mat.ptr(i, 0), mat.leadingDimension(),
                                     thread_id, stream_id);
}

// Scales the i[k]-th row of mat by val[k] for 0 <= k < i.size().
// In/Out: mat
// Preconditions: i.size() == val.size(), 0 <= i[k] < mat.nrRow() for 0 <= k < i.size().
template <typename ScalarType>
inline void scaleRows(Matrix<ScalarType, CPU>& mat, const Vector<int, CPU>& i,
                      const Vector<ScalarType, CPU>& val, int /*thread_id*/ = 0,
                      int /*stream_id*/ = 0) {
  assert(i.size() == val.size());

  for (int j = 0; j < mat.nrCols(); ++j)
    for (int ind = 0; ind < i.size(); ++ind)
      mat(i[ind], j) *= val[ind];
}
#ifdef DCA_HAVE_CUDA
template <typename ScalarType>
inline void scaleRows(Matrix<ScalarType, GPU>& mat, const Vector<int, GPU>& i,
                      const Vector<ScalarType, GPU>& val, int thread_id = 0, int stream_id = 0) {
  assert(i.size() == val.size());

  blas::scaleRows(mat.nrCols(), i.size(), i.ptr(), val.ptr(), mat.ptr(), mat.leadingDimension(),
                  thread_id, stream_id);
}
#endif  // DCA_HAVE_CUDA

// Swaps the j1-th column with the j2-th column of mat.
// In/Out: mat
// Preconditions: 0 <= j1 < mat.nrCols(), 0 <= j2 < mat_y.nrCols().
template <typename ScalarType, DeviceType device_name>
inline void swapCol(Matrix<ScalarType, device_name>& mat, int j1, int j2, int thread_id = 0,
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
template <typename ScalarType>
inline void swapCols(Matrix<ScalarType, GPU>& mat, const Vector<int, GPU>& j_1,
                     const Vector<int, GPU>& j_2, int thread_id = 0, int stream_id = 0) {
  assert(j_1.size() <= j_2.size());
  blas::swapCols(mat.nrRows(), j_1.size(), j_1.ptr(), j_2.ptr(), mat.ptr(), mat.leadingDimension(),
                 thread_id, stream_id);
}
#endif  // DCA_HAVE_CUDA

// Swaps the i1-th row with the i2-th row of mat.
// In/Out: mat
// Preconditions: 0 <= i1 < mat.nrRows(), 0 <= i2 < mat_y.nrRows().
template <typename ScalarType, DeviceType device_name>
inline void swapRow(Matrix<ScalarType, device_name>& mat, int i1, int i2, int thread_id = 0,
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
template <typename ScalarType>
inline void swapRows(Matrix<ScalarType, GPU>& mat, const Vector<int, GPU>& i_1,
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
template <typename ScalarType, DeviceType device_name>
inline void swapRowAndCol(Matrix<ScalarType, device_name>& mat, int i1, int i2, int thread_id = 0,
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
template <typename ScalarType>
void gemv(char transa, ScalarType alpha, const Matrix<ScalarType, CPU>& a,
          const Vector<ScalarType, CPU>& x, ScalarType beta, Vector<ScalarType, CPU>& y) {
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
template <typename ScalarType>
void gemv(char transa, const Matrix<ScalarType, CPU>& a, const Vector<ScalarType, CPU>& x,
          Vector<ScalarType, CPU>& y) {
  gemv<ScalarType>(transa, 1., a, x, 0., y);
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
template <typename ScalarType, DeviceType device_name, template <typename, DeviceType> class MatrixA,
          template <typename, DeviceType> class MatrixB, template <typename, DeviceType> class MatrixC>
void gemm(char transa, char transb, ScalarType alpha, const MatrixA<ScalarType, device_name>& a,
          const MatrixB<ScalarType, device_name>& b, ScalarType beta,
          MatrixC<ScalarType, device_name>& c, int thread_id = 0, int stream_id = 0) {
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
template <typename ScalarType, DeviceType device_name, template <typename, DeviceType> class MatrixA,
          template <typename, DeviceType> class MatrixB, template <typename, DeviceType> class MatrixC>
inline void gemm(const MatrixA<ScalarType, device_name>& a, const MatrixB<ScalarType, device_name>& b,
                 MatrixC<ScalarType, device_name>& c, int thread_id = 0, int stream_id = 0) {
  gemm<ScalarType, device_name>('N', 'N', 1., a, b, 0., c, thread_id, stream_id);
}

// Performs the matrix-matrix multiplication c <- alpha * a * b + beta * c,
// In/Out: c ('In' only if beta != 0)
// Preconditions: a.nrRows() == c.nrRows(), b.nrCols() == c.nrCols() and a.nrCols() == b.nrRows()
template <typename ScalarType, DeviceType device_name, template <typename, DeviceType> class MatrixA,
          template <typename, DeviceType> class MatrixB, template <typename, DeviceType> class MatrixC>
inline void gemm(ScalarType alpha, const MatrixA<ScalarType, device_name>& a,
                 const MatrixB<ScalarType, device_name>& b, ScalarType beta,
                 MatrixC<ScalarType, device_name>& c, int thread_id = 0, int stream_id = 0) {
  gemm<ScalarType, device_name>('N', 'N', alpha, a, b, beta, c, thread_id, stream_id);
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
template <typename ScalarType, DeviceType device_name>
inline void gemm(char transa, char transb, const Matrix<ScalarType, device_name>& a,
                 const Matrix<ScalarType, device_name>& b, Matrix<ScalarType, device_name>& c,
                 int thread_id = 0, int stream_id = 0) {
  gemm<ScalarType, device_name>(transa, transb, 1., a, b, 0., c, thread_id, stream_id);
}

// Performs the triangular solve b <- a^-1 * b,
// where a is a lower triangular matrix (uplo = 'L') or an upper triangular matrix (uplo = 'U'),
// with unit diagonal (diag = "U") or with general diagonal (diag = "N")
// In/Out: b
// Preconditions: a.nrRows() == a.nrCols() , a.nrCols() == b.nrRows()
template <typename ScalarType, DeviceType device_name>
void trsm(char uplo, char diag, const Matrix<ScalarType, device_name>& a,
          Matrix<ScalarType, device_name>& b, int thread_id = 0, int stream_id = 0) {
  assert(uplo == 'U' or uplo == 'L');
  assert(diag == 'U' or diag == 'N');
  assert(a.nrRows() == a.nrCols());
  assert(b.nrRows() == a.nrCols());

  blas::UseDevice<device_name>::trsm("L", &uplo, "N", &diag, b.nrRows(), b.nrCols(), ScalarType(1),
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
template <typename ScalarType>
void gemm(char transa, char transb, Matrix<ScalarType, CPU>& a,
          Matrix<std::complex<ScalarType>, CPU>& b, Matrix<std::complex<ScalarType>, CPU>& c) {
  Matrix<ScalarType, CPU> b_part(b.size());
  Matrix<ScalarType, CPU> c_re(c.size());
  Matrix<ScalarType, CPU> c_im(c.size());

  ScalarType sign = 1;
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

  gemm(transa, transb, sign, a, b_part, ScalarType(0), c_im);

  for (int j = 0; j < c.nrCols(); ++j)
    for (int i = 0; i < c.nrRows(); ++i)
      c(i, j) = std::complex<ScalarType>(c_re(i, j), c_im(i, j));
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
template <typename ScalarType>
static void gemm(char transa, char transb, Matrix<std::complex<ScalarType>, CPU>& a,
                 Matrix<ScalarType, CPU>& b, Matrix<std::complex<ScalarType>, CPU>& c) {
  Matrix<ScalarType, CPU> a_part(a.size());
  Matrix<ScalarType, CPU> c_re(c.size());
  Matrix<ScalarType, CPU> c_im(c.size());

  ScalarType sign = 1;
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

  gemm(transa, transb, sign, a_part, b, ScalarType(0), c_im);

  for (int j = 0; j < c.nrCols(); ++j)
    for (int i = 0; i < c.nrRows(); ++i)
      c(i, j) = std::complex<ScalarType>(c_re(i, j), c_im(i, j));
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
template <typename ScalarType>
void multiply(char transa, char transb, const std::array<Matrix<ScalarType, CPU>, 2>& a,
              const std::array<Matrix<ScalarType, CPU>, 2>& b,
              std::array<Matrix<ScalarType, CPU>, 2>& c,
              std::array<Matrix<ScalarType, CPU>, 5>& work) {
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

  const ScalarType signa = transa == 'C' ? transa = 'T', -1 : 1;
  const ScalarType signb = transb == 'C' ? transb = 'T', -1 : 1;

  for (int j = 0; j < a[0].nrCols(); ++j)
    for (int i = 0; i < a[0].nrRows(); ++i)
      a_sum(i, j) = a[0](i, j) + signa * a[1](i, j);
  for (int j = 0; j < b[0].nrCols(); ++j)
    for (int i = 0; i < b[0].nrRows(); ++i)
      b_sum(i, j) = b[0](i, j) + signb * b[1](i, j);

  gemm(transa, transb, a[0], b[0], work[0]);
  gemm(transa, transb, signa * signb, a[1], b[1], ScalarType(0), work[1]);
  gemm(transa, transb, a_sum, b_sum, work[2]);

  for (int j = 0; j < c[0].nrCols(); ++j)
    for (int i = 0; i < c[0].nrRows(); ++i) {
      c[0](i, j) = work[0](i, j) - work[1](i, j);
      c[1](i, j) = work[2](i, j) - work[0](i, j) - work[1](i, j);
    }
}

template <typename ScalarType>
void multiply(const std::array<Matrix<ScalarType, CPU>, 2>& a,
              const std::array<Matrix<ScalarType, CPU>, 2>& b,
              std::array<Matrix<ScalarType, CPU>, 2>& c,
              std::array<Matrix<ScalarType, CPU>, 5>& work) {
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
template <typename ScalarType, DeviceType device_name>
void multiply(char transa, char transb, const std::array<Matrix<ScalarType, device_name>, 2>& a,
              const Matrix<ScalarType, device_name>& b,
              std::array<Matrix<ScalarType, device_name>, 2>& c) {
  assert(transa == 'N' || transa == 'T' || transa == 'C');
  assert(transb == 'N' || transb == 'T');
  assert(a[0].size() == a[1].size());
  assert(c[0].size() == c[1].size());

  gemm(transa, transb, a[0], b, c[0]);
  const ScalarType sign = transa == 'C' ? transa = 'T', -1 : 1;
  gemm(transa, transb, sign, a[1], b, ScalarType(0), c[1]);
}

template <typename ScalarType, DeviceType device_name>
void multiply(const std::array<Matrix<ScalarType, device_name>, 2>& a,
              const Matrix<ScalarType, device_name>& b,
              std::array<Matrix<ScalarType, device_name>, 2>& c) {
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
template <typename ScalarType, DeviceType device_name>
void multiply(char transa, char transb, const Matrix<ScalarType, device_name>& a,
              const std::array<Matrix<ScalarType, device_name>, 2>& b,
              std::array<Matrix<ScalarType, device_name>, 2>& c) {
  assert(transa == 'N' || transa == 'T');
  assert(transb == 'N' || transb == 'T' || transb == 'C');
  assert(b[0].size() == b[1].size());
  assert(c[0].size() == c[1].size());

  gemm(transa, transb, a, b[0], c[0]);
  const ScalarType sign = transb == 'C' ? transb = 'T', -1 : 1;
  gemm(transa, transb, sign, a, b[1], ScalarType(0), c[1]);
}

template <typename ScalarType, DeviceType device_name>
void multiply(const Matrix<ScalarType, device_name>& a,
              const std::array<Matrix<ScalarType, device_name>, 2>& b,
              std::array<Matrix<ScalarType, device_name>, 2>& c) {
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
template <typename ScalarType, DeviceType device_name>
inline void multiplyDiagonalRight(const Matrix<ScalarType, device_name>& a,
                                  const Vector<ScalarType, device_name>& d,
                                  Matrix<ScalarType, device_name>& b, int thread_id = 0,
                                  int stream_id = 0) {
  lapack::UseDevice<device_name>::multiplyDiagonalRight(a.nrRows(), a.nrCols(), a.ptr(),
                                                        a.leadingDimension(), d.ptr(), 1, b.ptr(),
                                                        b.leadingDimension(), thread_id, stream_id);
}
template <typename ScalarType>
inline void multiplyDiagonalRight(const Matrix<ScalarType, GPU>& a, const Vector<ScalarType, CPU>& d,
                                  Matrix<ScalarType, GPU>& b, int thread_id = 0, int stream_id = 0) {
  Vector<ScalarType, GPU> d_gpu(d);
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
template <typename ScalarType>
void eigensolver(char jobvl, char jobvr, const Matrix<ScalarType, CPU>& a,
                 Vector<ScalarType, CPU>& lambda_re, Vector<ScalarType, CPU>& lambda_im,
                 Matrix<ScalarType, CPU>& vl, Matrix<ScalarType, CPU>& vr) {
  assert(a.is_square());

  Matrix<ScalarType, CPU> a_copy(a);
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
  dca::linalg::Vector<ScalarType, CPU> work(lwork);

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
template <typename ScalarType>
void eigensolver(char jobvl, char jobvr, const Matrix<std::complex<ScalarType>, CPU>& a,
                 Vector<std::complex<ScalarType>, CPU>& lambda,
                 Matrix<std::complex<ScalarType>, CPU>& vl,
                 Matrix<std::complex<ScalarType>, CPU>& vr) {
  assert(a.is_square());

  Matrix<std::complex<ScalarType>, CPU> a_copy(a);
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
  dca::linalg::Vector<std::complex<ScalarType>, CPU> work(lwork);
  dca::linalg::Vector<ScalarType, CPU> rwork(2 * a_copy.nrRows());

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
template <typename ScalarType>
void eigensolverSymmetric(char jobv, char uplo, const Matrix<ScalarType, CPU>& a,
                          Vector<ScalarType, CPU>& lambda, Matrix<ScalarType, CPU>& v) {
  assert(a.is_square());

  lambda.resizeNoCopy(a.nrRows());
  v = a;

  // Get optimal worksize.
  auto lwork = util::getEigensolverSymmetricWorkSize(jobv, uplo, v);
  dca::linalg::Vector<ScalarType, CPU> work(std::get<0>(lwork));
  dca::linalg::Vector<int, CPU> iwork(std::get<1>(lwork));

  lapack::syevd(&jobv, &uplo, v.nrRows(), v.ptr(), v.leadingDimension(), lambda.ptr(), work.ptr(),
                work.size(), iwork.ptr(), iwork.size());
}
// For real types Hermitian and symmetric is the same.
template <typename ScalarType>
inline void eigensolverHermitian(char jobv, char uplo, const Matrix<ScalarType, CPU>& a,
                                 Vector<ScalarType, CPU>& lambda, Matrix<ScalarType, CPU>& v) {
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
template <typename ScalarType>
void eigensolverHermitian(char jobv, char uplo, const Matrix<std::complex<ScalarType>, CPU>& a,
                          Vector<ScalarType, CPU>& lambda, Matrix<std::complex<ScalarType>, CPU>& v) {
  assert(a.is_square());

  lambda.resizeNoCopy(a.nrRows());
  v = a;

  // Get optimal worksize.
  auto lwork = util::getEigensolverHermitianWorkSize(jobv, uplo, v);
  dca::linalg::Vector<std::complex<ScalarType>, CPU> work(std::get<0>(lwork));
  dca::linalg::Vector<ScalarType, CPU> rwork(std::get<1>(lwork));
  dca::linalg::Vector<int, CPU> iwork(std::get<2>(lwork));

  lapack::heevd(&jobv, &uplo, v.nrRows(), v.ptr(), v.leadingDimension(), lambda.ptr(), work.ptr(),
                work.size(), rwork.ptr(), rwork.size(), iwork.ptr(), iwork.size());
}

template <typename ScalarType>
void eigensolverGreensFunctionMatrix(char jobv, char uplo,
                                     const Matrix<std::complex<ScalarType>, CPU>& a,
                                     Vector<ScalarType, CPU>& lambda,
                                     Matrix<std::complex<ScalarType>, CPU>& v) {
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
template <typename ScalarType>
void pseudoInverse(const Matrix<ScalarType, CPU>& a, Matrix<ScalarType, CPU>& a_inv,
                   double eps = 1.e-6) {
  int m = a.nrRows();
  int n = a.nrCols();
  a_inv.resizeNoCopy(std::make_pair(n, m));

  using RealType = decltype(std::real(*a.ptr()));

  if (m <= n) {
    // a_inv = a'*inv(a*a')
    // inv(a*a') = v*inv(lambda)*v', [lambda, v] = eig(a*a')

    Matrix<ScalarType, CPU> a_at("A_At", m);
    dca::linalg::matrixop::gemm('N', 'C', a, a, a_at);

    dca::linalg::Vector<RealType, CPU> lambda("Lambda", m);
    Matrix<ScalarType, CPU> v("V", m);

    eigensolverHermitian('V', 'U', a_at, lambda, v);
    Matrix<ScalarType, CPU> vt(v);

    for (int j = 0; j < m; j++) {
      ScalarType lambda_inv = 0;

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

    Matrix<ScalarType, CPU> at_a("at_a", n);
    dca::linalg::matrixop::gemm('C', 'N', a, a, at_a);

    dca::linalg::Vector<RealType, CPU> lambda("Lambda", n);
    Matrix<ScalarType, CPU> v("V", n);

    eigensolverHermitian('V', 'U', at_a, lambda, v);
    Matrix<ScalarType, CPU> vt(v);

    for (int j = 0; j < n; j++) {
      ScalarType lambda_inv = 0;

      if (lambda[j] > eps * lambda[n - 1])
        lambda_inv = 1. / lambda[j];

      scaleCol(v, j, lambda_inv);
    }

    gemm('N', 'C', v, vt, at_a);
    gemm('N', 'C', at_a, a, a_inv);
  }
}
}  // matrixop
}  // linalg
}  // dca

#endif  // DCA_LINALG_MATRIXOP_HPP
