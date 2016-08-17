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

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_FERMIONIC_ED_TYPE_DEFINITIONS_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_FERMIONIC_ED_TYPE_DEFINITIONS_H

#include <bitset>
#include <complex>

#include "comp_library/function_library/include_function_library.h"
#include "comp_library/linalg/linalg.hpp"

#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"

namespace DCA {
namespace ADVANCED_EXACT_DIAGONALIZATION {
// DCA::ADVANCED_EXACT_DIAGONALIZATION::

template <typename parameters_type>
struct advanced_ed_options {
public:
  const static std::size_t N = 8 * sizeof(std::size_t);

  using phi_type = std::bitset<N>;

  using profiler_t = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using int_type = int;
  using scalar_type = double;
  using complex_type = std::complex<scalar_type>;

  using vector_type = dca::linalg::Vector<scalar_type, dca::linalg::CPU>;
  using matrix_type = dca::linalg::Matrix<complex_type, dca::linalg::CPU>;
  using int_matrix_type = dca::linalg::Matrix<int, dca::linalg::CPU>;

  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index
  using nu_nu = dmn_variadic<nu, nu>;

  using r_DCA = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, CLUSTER,
                                     REAL_SPACE, BRILLOUIN_ZONE>>;
  using k_DCA = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, CLUSTER,
                                     MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  using b_dmn = b;
  using s_dmn = s;
  using r_dmn = r_DCA;
  using k_dmn = k_DCA;

  using nu_dmn = dmn_variadic<b_dmn, s_dmn>;
  using bs_dmn_type = dmn_variadic<b_dmn, s_dmn>;

  using nu_r_dmn_type = dmn_variadic<nu_dmn, r_dmn>;

  using b_s_r = dmn_variadic<b_dmn, s_dmn, r_dmn>;
  using bsr_dmn_type = dmn_variadic<b_dmn, s_dmn, r_dmn>;

  using b_s_k = dmn_variadic<b_dmn, s_dmn, k_dmn>;
  using bsk_dmn_type = dmn_variadic<b_dmn, s_dmn, k_dmn>;

  using nu_nu_r_dmn_type = dmn_variadic<nu_dmn, nu_dmn, r_dmn>;
  using nu_nu_k_dmn_type = dmn_variadic<nu_dmn, nu_dmn, k_dmn>;

public:
  static scalar_type get_epsilon() {
    return 1.e-3;
  }
};

}  // ADVANCED_EXACT_DIAGONALIZATION
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_FERMIONIC_ED_TYPE_DEFINITIONS_H
