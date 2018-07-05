// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Jérémie Bouquet (bouquetj@gmail.com).
//
//

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder_gpu.hpp"

#include "dca/linalg/make_constant_view.hpp"
#include "dca/linalg/matrix_view.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/kernels_interface.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/device_memory/global_memory_manager.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

using linalg::CPU;
using linalg::GPU;

DMatrixBuilder<linalg::GPU>::DMatrixBuilder(const G0Interpolation<GPU>& g0,
                                            const linalg::Matrix<int, linalg::CPU>& site_diff,
                                            const std::vector<int>& sbdm_step,
                                            const std::array<double, 3>& alpha)
    : BaseClass(g0.get_host_interpolation(), site_diff, sbdm_step, alpha), g0_ref_(g0) {
  ctint::GlobalMemoryManager::initializeCluster(*linalg::makeConstantView(site_diff), sbdm_step);
}

void DMatrixBuilder<GPU>::computeG0(Matrix& G0, const details::DeviceConfiguration& configuration,
                                    const int n_init, bool right_section, cudaStream_t stream) const {
  if (G0.nrRows() * G0.nrCols() == 0)
    return;
  details::buildG0Matrix(linalg::MatrixView<double, linalg::GPU>(G0), n_init, right_section,
                         configuration, g0_ref_, stream);
}

}  // ctint
}  // solver
}  // phys
}  // dca
