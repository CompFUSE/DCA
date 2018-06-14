// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
//
// Helper struct for ss_hybridization_walker_routines.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_STATIC_MATRIX_ROUTINES_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_STATIC_MATRIX_ROUTINES_HPP

#include <cassert>
#include <utility>

#include "dca/linalg/matrix.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

struct static_matrix_routines {
  template <typename scalar_type>
  static void cycle_column_forward(dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M) {
    assert(M.size().first == M.size().second);
    assert(M.capacity().first == M.capacity().second);

    std::pair<int, int> size = M.size();
    std::pair<int, int> capacity = M.capacity();

    scalar_type* tmp_column = new scalar_type[size.first];

    for (int l = 0; l < size.first; l++)
      tmp_column[l] = M(l, size.first - 1);

    if ((size.first - 1) > 0)
      memmove(&M(0, 1), &M(0, 0), sizeof(scalar_type) * capacity.first * (size.first - 1));

    for (int l = 0; l < size.first; l++)
      M(l, 0) = tmp_column[l];

    delete[] tmp_column;
  }

  template <typename scalar_type>
  static void cycle_column_backward(dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M) {
    assert(M.size().first == M.size().second);
    assert(M.capacity().first == M.capacity().second);

    std::pair<int, int> size = M.size();
    std::pair<int, int> capacity = M.capacity();

    scalar_type* tmp_column = new scalar_type[size.first];

    for (int l = 0; l < size.first; l++)
      tmp_column[l] = M(l, 0);

    if ((size.first - 1) > 0)
      memmove(&M(0, 0), &M(0, 1), sizeof(scalar_type) * capacity.first * (size.first - 1));

    for (int l = 0; l < size.first; l++)
      M(l, size.first - 1) = tmp_column[l];

    delete[] tmp_column;
  }

  template <typename scalar_type>
  static void cycle_row_forward(dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M) {
    assert(M.size().first == M.size().second);
    assert(M.capacity().first == M.capacity().second);

    std::pair<int, int> size = M.size();

    scalar_type* tmp_row = new scalar_type[size.first];

    for (int l = 0; l < size.first; l++)
      tmp_row[l] = M(size.first - 1, l);

    for (int l = 0; l < size.first; l++)
      for (int row_i = size.first - 1; row_i > 0; row_i--)
        M(row_i, l) = M(row_i - 1, l);

    for (int l = 0; l < size.first; l++)
      M(0, l) = tmp_row[l];

    delete[] tmp_row;
  }

  template <typename scalar_type>
  static void cycle_row_backward(dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M) {
    assert(M.size().first == M.size().second);
    assert(M.capacity().first == M.capacity().second);

    std::pair<int, int> size = M.size();

    scalar_type* tmp_row = new scalar_type[size.first];

    for (int l = 0; l < size.first; l++)
      tmp_row[l] = M(0, l);

    for (int l = 0; l < size.first; l++)
      for (int row_i = 0; row_i < size.first - 1; row_i++)
        M(row_i, l) = M(row_i + 1, l);

    for (int l = 0; l < size.first; l++)
      M(size.first - 1, l) = tmp_row[l];

    delete[] tmp_row;
  }
};

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_STATIC_MATRIX_ROUTINES_HPP
