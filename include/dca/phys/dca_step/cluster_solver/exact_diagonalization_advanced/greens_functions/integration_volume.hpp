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
// Helper class for TpGreensFunction.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_INTEGRATION_VOLUME_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_INTEGRATION_VOLUME_HPP

#include <vector>

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

template <typename scalar_type>
struct integration_volume_3D {
  bool is_converged;

  scalar_type x0_l, x0_m, x0_u, x1_l, x1_m, x1_u, x2_l, x2_m, x2_u;

  void do_recursion(std::vector<integration_volume_3D<scalar_type>>& /*vec*/) {
    std::vector<integration_volume_3D<scalar_type>> vols(8);

    int index = 0;
    for (int l0 = 0; l0 < 2; l0++) {
      for (int l1 = 0; l1 < 2; l1++) {
        for (int l2 = 0; l2 < 2; l2++) {
          integration_volume_3D<scalar_type>& volume = vols[index];

          volume.x0_l = l0 == 0 ? x0_l : x0_m;
          volume.x0_u = l0 == 0 ? x0_m : x0_u;

          volume.x1_l = l1 == 0 ? x1_l : x1_m;
          volume.x1_u = l1 == 0 ? x1_m : x1_u;

          volume.x2_l = l2 == 0 ? x2_l : x2_m;
          volume.x2_u = l2 == 0 ? x2_m : x2_u;

          index += 1;
        }
      }
    }

    for (int l0 = 0; l0 < 8; l0++) {
      integration_volume_3D<scalar_type>& volume = vols[l0];

      volume.x0_m = (x0_l + x0_u) / 2.;
      volume.x1_m = (x1_l + x1_u) / 2.;
      volume.x2_m = (x2_l + x2_u) / 2.;
    }
  }
};

template <typename scalar_type>
struct integration_volume_4D {
  bool is_converged;

  scalar_type x0_l, x0_m, x0_u, x1_l, x1_m, x1_u, x2_l, x2_m, x2_u, x3_l, x3_m, x3_u;

  void do_recursion(std::vector<integration_volume_4D<scalar_type>>& /*vec*/) {}
};

}  // ed
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_INTEGRATION_VOLUME_HPP
