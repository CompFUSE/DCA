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
// This file provides type definitions and utility functions for the ED solver.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_OPTIONS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_OPTIONS_HPP

#include <bitset>
#include <complex>

#include "dca/function/domains.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

#include "comp_library/linalg/linalg.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

template <typename parameters_type>
struct Options {
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

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index
  using nu_nu = func::dmn_variadic<nu, nu>;

  using r_DCA =
      func::dmn_0<domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                                          domains::REAL_SPACE, domains::BRILLOUIN_ZONE>>;
  using k_DCA =
      func::dmn_0<domains::cluster_domain<double, parameters_type::lattice_type::DIMENSION, domains::CLUSTER,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;

  using b_dmn = b;
  using s_dmn = s;
  using r_dmn = r_DCA;
  using k_dmn = k_DCA;

  using nu_dmn = func::dmn_variadic<b_dmn, s_dmn>;
  using bs_dmn_type = func::dmn_variadic<b_dmn, s_dmn>;

  using nu_r_dmn_type = func::dmn_variadic<nu_dmn, r_dmn>;

  using b_s_r = func::dmn_variadic<b_dmn, s_dmn, r_dmn>;
  using bsr_dmn_type = func::dmn_variadic<b_dmn, s_dmn, r_dmn>;

  using b_s_k = func::dmn_variadic<b_dmn, s_dmn, k_dmn>;
  using bsk_dmn_type = func::dmn_variadic<b_dmn, s_dmn, k_dmn>;

  using nu_nu_r_dmn_type = func::dmn_variadic<nu_dmn, nu_dmn, r_dmn>;
  using nu_nu_k_dmn_type = func::dmn_variadic<nu_dmn, nu_dmn, k_dmn>;

  static scalar_type get_epsilon() {
    return 1.e-3;
  }
};

}  // ed
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_OPTIONS_HPP
