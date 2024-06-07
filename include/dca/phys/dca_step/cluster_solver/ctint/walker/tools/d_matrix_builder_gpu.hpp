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

#include "dca/config/haves_defines.hpp"

#ifndef DCA_HAVE_GPU
#error "This file needs GPU support."
#endif

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"

#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration.hpp"
#include "dca/linalg/util/gpu_stream.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/g0_interpolation_gpu.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <typename Scalar>
class DMatrixBuilder<linalg::GPU, Scalar> final : public DMatrixBuilder<linalg::CPU, Scalar> {
private:
  using Matrix = linalg::Matrix<Scalar, linalg::GPU>;
  using MatrixPair = std::array<Matrix, 2>;
  using BaseClass = DMatrixBuilder<linalg::CPU, Scalar>;
  using GpuStream = dca::linalg::util::GpuStream;
public:
  template <class RDmn>
  DMatrixBuilder(const G0Interpolation<linalg::GPU, Scalar>& g0, int nb, const RDmn& /*rdmn*/);

  DMatrixBuilder(const G0Interpolation<linalg::GPU, Scalar>& g0,
                 const linalg::Matrix<int, linalg::CPU>& site_diff,
                 const linalg::Matrix<int, linalg::CPU>& site_add, int nb, int r0);

  const auto& getG0() const {
    return g0_ref_;
  }

  // See DMatrixBuilder<linalg::CPU, Scalar>::computeG0.
  // Out: G0. Device matrix
  void computeG0(Matrix& G0, const details::DeviceConfiguration& configuration, int n_init,
                 bool right_section, const GpuStream& stream) const;

private:
  const G0Interpolation<linalg::GPU, Scalar>& g0_ref_;
};

template <typename Scalar>
template <class RDmn>
DMatrixBuilder<linalg::GPU, Scalar>::DMatrixBuilder(const G0Interpolation<linalg::GPU, Scalar>& g0,
                                                  int nb, const RDmn& /*rdmn*/)
    : DMatrixBuilder(g0, RDmn::parameter_type::get_add_matrix(),
                     RDmn::parameter_type::get_subtract_matrix(), nb,
                     RDmn::parameter_type::origin_index()) {}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_D_MATRIX_BUILDER_GPU_HPP
