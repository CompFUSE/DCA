// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//          Jérémie Bouquet (bouquetj@gmail.com)
//
// Builds G0 matrices used by the GPU walker.

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
class DMatrixBuilder<linalg::GPU> final : public DMatrixBuilder<linalg::CPU> {
private:
  using Matrix = linalg::Matrix<double, linalg::GPU>;
  using MatrixPair = std::array<linalg::Matrix<double, linalg::GPU>, 2>;
  using BaseClass = DMatrixBuilder<linalg::CPU>;

public:
  template <class RDmn>
  DMatrixBuilder(const G0Interpolation<linalg::GPU>& g0, int nb, const RDmn& /*rdmn*/);

  DMatrixBuilder(const G0Interpolation<linalg::GPU>& g0,
                 const linalg::Matrix<int, linalg::CPU>& site_diff,
                 const linalg::Matrix<int, linalg::CPU>& site_add, int nb, int r0);


  const G0Interpolation<linalg::GPU>& getG0() const {
    return g0_ref_;
  }

  void computeG0(Matrix& G0, const details::DeviceConfiguration& configuration, int n_init,
                 bool right_section, cudaStream_t stream) const override;

private:
  const G0Interpolation<linalg::GPU>& g0_ref_;
};

template <class RDmn>
DMatrixBuilder<linalg::GPU>::DMatrixBuilder(const G0Interpolation<linalg::GPU>& g0, int nb,
                                            const RDmn& /*rdmn*/)
    : DMatrixBuilder(g0, RDmn::parameter_type::get_add_matrix(),
                     RDmn::parameter_type::get_subtract_matrix(), nb,
                     RDmn::parameter_type::origin_index()) {}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_GPU_HPP
