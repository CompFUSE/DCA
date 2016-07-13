// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Description

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_OVERLAP_MATRIX_ELEMENT_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_OVERLAP_MATRIX_ELEMENT_H

namespace DCA {
namespace ADVANCED_EXACT_DIAGONALIZATION {
// DCA::ADVANCED_EXACT_DIAGONALIZATION::

template <typename parameter_type, typename ed_options>
struct sparse_element {
  typedef typename ed_options::scalar_type scalar_type;
  typedef typename ed_options::complex_type complex_type;

  int i, j;
  complex_type value;
};

template <typename parameter_type, typename ed_options>
bool operator<(const sparse_element<parameter_type, ed_options /*b_dmn, s_dmn, r_dmn*/>& el_1,
               const sparse_element<parameter_type, ed_options /*b_dmn, s_dmn, r_dmn*/>& el_2) {
  if (el_1.j < el_2.j)
    return true;

  else if (el_1.j == el_2.j && el_1.i < el_2.i)
    return true;

  else
    return false;
}

}  // ADVANCED_EXACT_DIAGONALIZATION
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_OVERLAP_MATRIX_ELEMENT_H
