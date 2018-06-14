// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides a helper class for fermionic_overlap_matrices to store the elements of the
// sparse matrices.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_OVERLAP_MATRIX_ELEMENT_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_OVERLAP_MATRIX_ELEMENT_HPP

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

template <typename parameter_type, typename ed_options>
struct OverlapMatrixElement {
  typedef typename ed_options::scalar_type scalar_type;
  typedef typename ed_options::complex_type complex_type;

  int i, j;
  complex_type value;
};

template <typename parameter_type, typename ed_options>
bool operator<(const OverlapMatrixElement<parameter_type, ed_options>& el_1,
               const OverlapMatrixElement<parameter_type, ed_options>& el_2) {
  if (el_1.j < el_2.j)
    return true;

  else if (el_1.j == el_2.j && el_1.i < el_2.i)
    return true;

  else
    return false;
}

}  // ed
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_OVERLAP_MATRIX_ELEMENT_HPP
