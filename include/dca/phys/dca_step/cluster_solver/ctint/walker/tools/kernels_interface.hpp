// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Kernel interface for GPU interpolation.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_INTERFACE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_INTERFACE_HPP
#ifdef DCA_HAVE_GPU

#include "dca/platform/dca_gpu.h"
#include "dca/linalg/util/gpu_stream.hpp"

#include "dca/linalg/matrix_view.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/device_interpolation_data.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::
  
template <typename Scalar, typename SignType>
void buildG0Matrix(linalg::MatrixView<Scalar, linalg::GPU> G0, const int n_init,
                   const bool right_section, DeviceConfiguration config,
                   DeviceInterpolationData<Scalar, SignType> g0_interp, const dca::linalg::util::GpuStream& stream);
extern  template void buildG0Matrix(linalg::MatrixView<float, linalg::GPU>, const int, const bool,
                            DeviceConfiguration, DeviceInterpolationData<float, signed char>, const dca::linalg::util::GpuStream&);
extern  template void buildG0Matrix(linalg::MatrixView<double, linalg::GPU>, const int, const bool,
                            DeviceConfiguration, DeviceInterpolationData<double, std::int8_t>, const dca::linalg::util::GpuStream&);
extern  template void buildG0Matrix(linalg::MatrixView<std::complex<float>, linalg::GPU>, const int,
                            const bool, DeviceConfiguration,
                            DeviceInterpolationData<std::complex<float>, std::complex<float>>, const dca::linalg::util::GpuStream&);
extern  template void buildG0Matrix(linalg::MatrixView<std::complex<double>, linalg::GPU>, const int,
                            const bool, DeviceConfiguration,
                            DeviceInterpolationData<std::complex<double>, std::complex<double>>, const dca::linalg::util::GpuStream&);

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_HAVE_GPU
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_KERNELS_INTERFACE_HPP
