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
// - insertCol, insertRow (for CPU matrices only)
// - removeCol, removeRow, removeRowAndCol
// - scaleCol, scaleRow, scaleRows
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
