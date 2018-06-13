// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Jérémie Bouquet (bouquetj@gmail.com)
//
//

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_GPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_GPU_HPP

#ifndef DCA_HAVE_CUDA
#error "This file needs GPU support."
#endif

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"

#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation_gpu.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <>
class DMatrixBuilder<linalg::GPU> : public DMatrixBuilder<linalg::CPU> {
private:
  using Matrix = linalg::Matrix<double, linalg::GPU>;
  using MatrixPair = std::array<linalg::Matrix<double, linalg::GPU>, 2>;
  using BaseClass = DMatrixBuilder<linalg::CPU>;

public:
  DMatrixBuilder(const G0Interpolation<linalg::GPU>& g0,
                 const linalg::Matrix<int, linalg::CPU>& site_diff,
                 const std::vector<int>& sbdm_step, const std::array<double, 3>& alphas);

  const G0Interpolation<linalg::GPU>& getG0() const {
    return g0_ref_;
  }

  void computeG0(Matrix& G0, const details::DeviceConfiguration& configuration, int n_init,
                 bool right_section, cudaStream_t stream) const;

private:
  const G0Interpolation<linalg::GPU>& g0_ref_;
};

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_GPU_HPP
