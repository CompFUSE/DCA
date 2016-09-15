// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the matrix interface for the following matrix operations:
// - copyCol, copyRow, copyCols, copyRows
// - difference
// - insertCol, insertRow (for CPU matrices only)
// - removeCol, removeRow, removeRowAndCol
// - scaleCol, scaleRow, scaleRows
// - swapCol, swapRow, swapRowAndCol
// - swapCols, swapRows (for GPU matrices only)
// - gemm
// - trsm

#ifndef DCA_LINALG_MATRIXOP_HPP
#define DCA_LINALG_MATRIXOP_HPP

#include <cassert>
#include <cstring>

#include "dca/linalg/blas/use_device.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

#include "comp_library/linalg/src/linalg_operations/GPUfunc.h"

namespace dca {
namespace linalg {
namespace matrixop {
// dca::linalg::matrixop::

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
template <typename ScalarType>
inline void copyCols(const Matrix<ScalarType, GPU>& mat_x, const Vector<int, GPU>& j_x,
                     Matrix<ScalarType, GPU>& mat_y, const Vector<int, GPU>& j_y, int thread_id = 0,
                     int stream_id = 0) {
  assert(j_x.size() <= j_y.size());
  assert(mat_x.nrRows() == mat_y.nrRows());

  GPU_KERNEL::many_column_copies(mat_x.nrRows(), j_x.size(), j_x.ptr(), mat_x.ptr(),
                                 mat_x.leadingDimension(), j_y.ptr(), mat_y.ptr(),
                                 mat_y.leadingDimension(), thread_id, stream_id);
}

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
template <typename ScalarType>
inline void copyRows(const Matrix<ScalarType, GPU>& mat_x, const Vector<int, GPU>& i_x,
                     Matrix<ScalarType, GPU>& mat_y, const Vector<int, GPU>& i_y, int thread_id = 0,
                     int stream_id = 0) {
  assert(i_x.size() <= i_y.size());
  assert(mat_x.nrCols() == mat_y.nrCols());

  GPU_KERNEL::many_row_copies(mat_x.nrCols(), i_x.size(), i_x.ptr(), mat_x.ptr(),
                              mat_x.leadingDimension(), i_y.ptr(), mat_y.ptr(),
                              mat_y.leadingDimension(), thread_id, stream_id);
}

// Returns the difference of two matrices in terms of max_i,j(|a(i, j) - b(i, j)|).
// If the difference is larger than the threshold a std::logic_error exception is thrown,
// and if DNDEBUG is not defined each difference which exceeds the threshold is printed.
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
#ifndef DNDEBUG
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
#endif  // DNDEBUG

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
            sizeof(ScalarType) * (mat.nrCols() - j) * mat.leadingDimension());

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
      memmove(mat.ptr(i + 1, j), mat.ptr(i, j), sizeof(ScalarType) * (mat.nrRows() - i));

  for (int j = 0; j < mat.nrCols(); ++j)
    mat(i, j) = 0;
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
template <typename ScalarType>
void removeCol(Matrix<ScalarType, GPU>& mat, int j) {
  assert(j >= 0 && j < mat.nrCols());

  if (mat.nrRows() > 0 && j < mat.nrCols() - 1)
    LIN_ALG::MEMORY_MANAGEMENT<GPU>::remove_first_col(mat.nrRows(), mat.nrCols() - j, mat.ptr(0, j),
                                                      mat.leadingDimension());

  mat.resize(std::make_pair(mat.nrRows(), mat.nrCols() - 1));
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
template <typename ScalarType>
void removeRow(Matrix<ScalarType, GPU>& mat, int i) {
  assert(i >= 0 && i < mat.nrRows());

  if (mat.nrCols() > 0 && i < mat.nrRows() - 1)
    LIN_ALG::MEMORY_MANAGEMENT<GPU>::remove_first_row(mat.nrRows() - i, mat.nrCols(), mat.ptr(i, 0),
                                                      mat.leadingDimension());

  mat.resize(std::make_pair(mat.nrRows() - 1, mat.nrCols()));
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
template <typename ScalarType>
inline void scaleRows(Matrix<ScalarType, GPU>& mat, const Vector<int, GPU>& i,
                      const Vector<ScalarType, GPU>& val, int thread_id = 0, int stream_id = 0) {
  assert(i.size() == val.size());

  GPU_KERNEL::scale_many_rows(mat.nrCols(), i.size(), i.ptr(), val.ptr(), mat.ptr(),
                              mat.leadingDimension(), thread_id, stream_id);
}

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
template <typename ScalarType>
inline void swapCols(Matrix<ScalarType, GPU>& mat, const Vector<int, GPU>& j_1,
                     const Vector<int, GPU>& j_2, int thread_id = 0, int stream_id = 0) {
  assert(j_1.size() <= j_2.size());
  LIN_ALG::GPU_KERNEL::swap_many_cols(mat.nrRows(), mat.nrCols(), mat.ptr(), mat.leadingDimension(),
                                      j_1.size(), j_1.ptr(), j_2.ptr(), thread_id, stream_id);
}

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
template <typename ScalarType>
inline void swapRows(Matrix<ScalarType, GPU>& mat, const Vector<int, GPU>& i_1,
                     const Vector<int, GPU>& i_2, int thread_id = 0, int stream_id = 0) {
  assert(i_1.size() == i_2.size());
  LIN_ALG::GPU_KERNEL::swap_many_rows(mat.nrRows(), mat.nrCols(), mat.ptr(), mat.leadingDimension(),
                                      i_1.size(), i_1.ptr(), i_2.ptr(), thread_id, stream_id);
}

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

// Performs the matrix-matrix multiplication c <- alpha * op(a) * op(b) + beta * c,
// where op(X) = X if transX == 'N', op(X) = transposed(X) if transX == 'T', and
// op(X) == conjugate_transposed(X) if transX == 'C' (X = a, b).
// In/Out: c ('In' only if beta != 0)
// Preconditions: transa and transb should be one of the following: 'N', 'T' for real ScalarType
//                or  'N', 'T', 'C' for complex ScalarType
//                a.nrRows() == c.nrRows() if transa == 'N', a.nrCols() == c.nrRows() otherwise,
//                b.nrCols() == c.nrCols() if transb == 'N', b.nrRows() == c.nrCols() otherwise,
//                ka == kb, where ka = a.nrCols() if transa == 'N', ka = a.nrRows() otherwise and
//                          kb = b.nrRows() if transb == 'N', kb = b.nrCols() otherwise.
template <typename ScalarType, DeviceType device_name>
void gemm(char transa, char transb, ScalarType alpha, const matrix<ScalarType, device_name>& a,
          const matrix<ScalarType, device_name>& b, ScalarType beta,
          matrix<ScalarType, device_name>& c, int thread_id = 0, int stream_id = 0) {
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
template <typename ScalarType, DeviceType device_name>
inline void gemm(const matrix<ScalarType, device_name>& a, const matrix<ScalarType, device_name>& b,
                 matrix<ScalarType, device_name>& c, int thread_id = 0, int stream_id = 0) {
  gemm<ScalarType, device_name>('N', 'N', 1., a, b, 0., c, thread_id, stream_id);
}

// Performs the matrix-matrix multiplication c <- alpha * a * b + beta * c,
// In/Out: c ('In' only if beta != 0)
// Preconditions: a.nrRows() == c.nrRows(), b.nrCols() == c.nrCols() and a.nrCols() == b.nrRows()
template <typename ScalarType, DeviceType device_name>
inline void gemm(ScalarType alpha, const matrix<ScalarType, device_name>& a,
                 const matrix<ScalarType, device_name>& b, ScalarType beta,
                 matrix<ScalarType, device_name>& c, int thread_id = 0, int stream_id = 0) {
  gemm<ScalarType, device_name>('N', 'N', alpha, a, b, beta, c, thread_id, stream_id);
}

// Performs the matrix-matrix multiplication c <- op(a) * op(b),
// where op(X) = X if transX == 'N', op(X) = transposed(X) if transX == 'T', and
// op(X) == conjugate_transposed(X) if transX == 'C' (X = a, b).
// Out: c
// Preconditions: transa and transb should be one of the following: 'N', 'T' for real ScalarType
//                or  'N', 'T', 'C' for complex ScalarType
//                a.nrRows() == c.nrRows() if transa == 'N', a.nrCols() == c.nrRows() otherwise,
//                b.nrCols() == c.nrCols() if transb == 'N', b.nrRows() == c.nrCols() otherwise,
//                ka == kb, where ka = a.nrCols() if transa == 'N', ka = a.nrRows() otherwise and
//                          kb = b.nrRows() if transb == 'N', kb = b.nrCols() otherwise.
template <typename ScalarType, DeviceType device_name>
inline void gemm(char transa, char transb, const matrix<ScalarType, device_name>& a,
                 const matrix<ScalarType, device_name>& b, matrix<ScalarType, device_name>& c,
                 int thread_id = 0, int stream_id = 0) {
  gemm<ScalarType, device_name>(transa, transb, 1., a, b, 0., c, thread_id, stream_id);
}

// Performs the triangular solve b <- a^-1 * b,
// where a is a lower triangular matrix (uplo = 'L') or an upper triangular matrix (uplo = 'U'),
// with unit diagonal (diag = "U") or with general diagonal (diag = "N")
// In/Out: b
// Preconditions: a.nrRows() == a.nrCols() , a.nrCols() == b.nrRows()
template <typename ScalarType, DeviceType device_name>
void trsm(char uplo, char diag, const matrix<ScalarType, device_name>& a,
          matrix<ScalarType, device_name>& b, int thread_id = 0, int stream_id = 0) {
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
  matrix<ScalarType, CPU> b_part(b.size());
  matrix<ScalarType, CPU> c_re(c.size());
  matrix<ScalarType, CPU> c_im(c.size());

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
  matrix<ScalarType, CPU> a_part(a.size());
  matrix<ScalarType, CPU> c_re(c.size());
  matrix<ScalarType, CPU> c_im(c.size());

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

}  // matrixop
}  // linalg
}  // dca

#endif  // DCA_LINALG_MATRIXOP_HPP
