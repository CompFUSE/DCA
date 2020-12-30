// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class organizes the interpolation G0(tau) for tau in [0, beta]
// specialization for CPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_WALKER_METHODS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_WALKER_METHODS_HPP

#include "dca/linalg/matrix_view.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::

enum MoveType { VERTEX_INSERTION, VERTEX_REMOVAL };

inline double computeAcceptanceProbability(const int order_change, const double det_ratio,
                                           const double total_interaction, const double beta,
                                           const double w, const int combinatorial_factor,
                                           MoveType type) {
  double strength_factor;
  switch (order_change) {
    case 1:
      strength_factor = beta * (w > 0 ? -total_interaction : total_interaction);
      break;
    case 2:
      strength_factor = beta * beta * std::abs(w) * total_interaction;
      break;
    default:
      throw(std::logic_error("General case not supported"));
  }
  return type == VERTEX_INSERTION ? det_ratio * strength_factor / combinatorial_factor
                                  : det_ratio * combinatorial_factor / strength_factor;
}

template <class MatrixA>
inline void smallInverse(MatrixA& m) {
  assert(m.is_square());
  switch (m.nrCols()) {
    case 0:
      return;
    case 1:
      m(0, 0) = 1. / m(0, 0);
      break;
    case 2: {
      const double det = m(1, 1) * m(0, 0) - m(1, 0) * m(0, 1);
      std::swap(m(0, 0), m(1, 1));
      m(0, 0) /= det;
      m(1, 0) /= -det;
      m(0, 1) /= -det;
      m(1, 1) /= det;
    } break;
    default:
      linalg::matrixop::inverse(m);
  }
}

template <class MatrixA, class MatrixB>
inline void smallInverse(const MatrixA& in, MatrixB& out) {
  assert(in.size() == out.size());
  switch (in.nrCols()) {
    case 1:
      out(0, 0) = 1. / in(0, 0);
      break;
    case 2: {
      const double det = in(1, 1) * in(0, 0) - in(1, 0) * in(0, 1);
      out(0, 0) = in(1, 1) / det;
      out(0, 1) = -in(0, 1) / det;
      out(1, 0) = -in(1, 0) / det;
      out(1, 1) = in(0, 0) / det;
    } break;
    default:
      out = in;
      linalg::matrixop::inverse(out);
  }
}

template <class MatrixA, class MatrixB, typename Real>
inline void smallInverse(const MatrixA& in, MatrixB& out, const Real det,
                         linalg::Vector<int, linalg::CPU>& ipiv,
                         linalg::Vector<Real, linalg::CPU>& work) {
  static_assert(std::is_same<std::remove_cv_t<typename MatrixA::ValueType>, Real>::value,
                "Scalar type MatrixA mismatch.");
  static_assert(std::is_same<typename MatrixB::ValueType, Real>::value,
                "Scalar type MatrixB mismatch.");

  assert(in.size() == out.size());
  switch (in.nrCols()) {
    case 1:
      out(0, 0) = 1. / in(0, 0);
      break;
    case 2:
      out(0, 0) = in(1, 1) / det;
      out(0, 1) = -in(0, 1) / det;
      out(1, 0) = -in(1, 0) / det;
      out(1, 1) = in(0, 0) / det;
      break;
    default:
      // TODO improve
      out = in;
      linalg::matrixop::inverse(out, ipiv, work);
  }
}

template <class Matrix, typename Real>
inline void smallInverse(Matrix& m, const Real det, linalg::Vector<int, linalg::CPU>& ipiv,
                         linalg::Vector<Real, linalg::CPU>& work) {
  static_assert(std::is_same<typename Matrix::ValueType, Real>::value, "Matrix type mismatch.");
  assert(m.size() == m.size());

  switch (m.nrCols()) {
    case 1:
      m(0, 0) = 1. / m(0, 0);
      break;
    case 2:
      assert(det);
      {
        const Real tmp = m(0, 0);
        m(0, 0) = m(1, 1) / det;
        m(1, 0) /= -det;
        m(0, 1) /= -det;
        m(1, 1) = tmp / det;
      }
      break;
    default:
      linalg::matrixop::inverse(m, ipiv, work);
  }
}

template <class MatrixType>
inline double smallDeterminant(const MatrixType& M) {
  static_assert(MatrixType::device == linalg::CPU, "GPU small inverse is not defined.");
  assert(M.is_square());
  switch (M.nrCols()) {
    case 0:
      return 1;
    case 1:
      return M(0, 0);
    case 2:
      return M(0, 0) * M(1, 1) - M(1, 0) * M(0, 1);
    case 3:
      return M(0, 0) * M(1, 1) * M(2, 2) + M(0, 1) * M(1, 2) * M(2, 0) + M(0, 2) * M(1, 0) * M(2, 1) -
             M(0, 2) * M(1, 1) * M(2, 0) - M(2, 1) * M(1, 2) * M(0, 0) - M(2, 2) * M(1, 0) * M(0, 1);
    case 4:
      return M(0, 0) * (M(1, 1) * M(2, 2) * M(3, 3) + M(1, 2) * M(2, 3) * M(3, 1) +
                        M(1, 3) * M(2, 1) * M(3, 2) - M(1, 3) * M(2, 2) * M(3, 1) -
                        M(3, 2) * M(2, 3) * M(1, 1) - M(3, 3) * M(2, 1) * M(1, 2)) -
             M(0, 1) * (M(1, 0) * M(2, 2) * M(3, 3) + M(1, 2) * M(2, 3) * M(3, 0) +
                        M(1, 3) * M(2, 0) * M(3, 2) - M(1, 3) * M(2, 2) * M(3, 0) -
                        M(3, 2) * M(2, 3) * M(1, 0) - M(3, 3) * M(2, 0) * M(1, 2)) +
             M(0, 2) * (M(1, 0) * M(2, 1) * M(3, 3) + M(1, 1) * M(2, 3) * M(3, 0) +
                        M(1, 3) * M(2, 0) * M(3, 1) - M(1, 3) * M(2, 1) * M(3, 0) -
                        M(3, 1) * M(2, 3) * M(1, 0) - M(3, 3) * M(2, 0) * M(1, 1)) -
             M(0, 3) * (M(1, 0) * M(2, 1) * M(3, 2) + M(1, 1) * M(2, 2) * M(3, 0) +
                        M(1, 2) * M(2, 0) * M(3, 1) - M(1, 2) * M(2, 1) * M(3, 0) -
                        M(3, 1) * M(2, 2) * M(1, 0) - M(3, 2) * M(2, 0) * M(1, 1));
    default:
      throw(std::logic_error("General case not supported."));
  }
}

template <class Matrix, class Vector>
inline double separateIndexDeterminant(const Matrix& M, const Vector& indices) {
  switch (indices.size()) {
    case 1:
      return M(indices[0], indices[0]);
    case 2:
      return M(indices[0], indices[0]) * M(indices[1], indices[1]) -
             M(indices[1], indices[0]) * M(indices[0], indices[1]);
    case 3:
      return M(indices[0], indices[0]) * M(indices[1], indices[1]) * M(indices[2], indices[2]) +
             M(indices[0], indices[1]) * M(indices[1], indices[2]) * M(indices[2], indices[0]) +
             M(indices[0], indices[2]) * M(indices[1], indices[0]) * M(indices[2], indices[1]) -
             M(indices[0], indices[2]) * M(indices[1], indices[1]) * M(indices[2], indices[0]) -
             M(indices[2], indices[1]) * M(indices[1], indices[2]) * M(indices[0], indices[0]) -
             M(indices[2], indices[2]) * M(indices[1], indices[0]) * M(indices[0], indices[1]);
    case 4:
      return M(indices[0], indices[0]) *
                 (M(indices[1], indices[1]) * M(indices[2], indices[2]) * M(indices[3], indices[3]) +
                  M(indices[1], indices[2]) * M(indices[2], indices[3]) * M(indices[3], indices[1]) +
                  M(indices[1], indices[3]) * M(indices[2], indices[1]) * M(indices[3], indices[2]) -
                  M(indices[1], indices[3]) * M(indices[2], indices[2]) * M(indices[3], indices[1]) -
                  M(indices[3], indices[2]) * M(indices[2], indices[3]) * M(indices[1], indices[1]) -
                  M(indices[3], indices[3]) * M(indices[2], indices[1]) * M(indices[1], indices[2])) -
             M(indices[0], indices[1]) *
                 (M(indices[1], indices[0]) * M(indices[2], indices[2]) * M(indices[3], indices[3]) +
                  M(indices[1], indices[2]) * M(indices[2], indices[3]) * M(indices[3], indices[0]) +
                  M(indices[1], indices[3]) * M(indices[2], indices[0]) * M(indices[3], indices[2]) -
                  M(indices[1], indices[3]) * M(indices[2], indices[2]) * M(indices[3], indices[0]) -
                  M(indices[3], indices[2]) * M(indices[2], indices[3]) * M(indices[1], indices[0]) -
                  M(indices[3], indices[3]) * M(indices[2], indices[0]) * M(indices[1], indices[2])) +
             M(indices[0], indices[2]) *
                 (M(indices[1], indices[0]) * M(indices[2], indices[1]) * M(indices[3], indices[3]) +
                  M(indices[1], indices[1]) * M(indices[2], indices[3]) * M(indices[3], indices[0]) +
                  M(indices[1], indices[3]) * M(indices[2], indices[0]) * M(indices[3], indices[1]) -
                  M(indices[1], indices[3]) * M(indices[2], indices[1]) * M(indices[3], indices[0]) -
                  M(indices[3], indices[1]) * M(indices[2], indices[3]) * M(indices[1], indices[0]) -
                  M(indices[3], indices[3]) * M(indices[2], indices[0]) * M(indices[1], indices[1])) -
             M(indices[0], indices[3]) *
                 (M(indices[1], indices[0]) * M(indices[2], indices[1]) * M(indices[3], indices[2]) +
                  M(indices[1], indices[1]) * M(indices[2], indices[2]) * M(indices[3], indices[0]) +
                  M(indices[1], indices[2]) * M(indices[2], indices[0]) * M(indices[3], indices[1]) -
                  M(indices[1], indices[2]) * M(indices[2], indices[1]) * M(indices[3], indices[0]) -
                  M(indices[3], indices[1]) * M(indices[2], indices[2]) * M(indices[1], indices[0]) -
                  M(indices[3], indices[2]) * M(indices[2], indices[0]) * M(indices[1], indices[1]));
    default:
      throw(std::logic_error("General case not supported."));
  }
}

template <typename Scalar>
inline void removeIndices(linalg::Matrix<Scalar, linalg::CPU>& Q,
                          linalg::Matrix<Scalar, linalg::CPU>& R,
                          const std::vector<unsigned short>& indices) {
  if (indices.size() == 0)
    return;

  int destination = Q.nrRows() - 1;
  for (int vertex : indices) {
    // TODO single line
    linalg::matrixop::swapRow(Q, vertex, destination);

    linalg::matrixop::swapCol(R, vertex, destination);

    --destination;
  }

  Q.resize(std::make_pair(Q.nrRows() - indices.size(), Q.nrCols()));
  R.resize(std::make_pair(R.nrRows(), R.nrCols() - indices.size()));
}

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_WALKER_METHODS_HPP
