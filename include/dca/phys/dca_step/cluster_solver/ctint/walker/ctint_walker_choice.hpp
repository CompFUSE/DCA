// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Selects the appropriate ctint walker implementation.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CHOICE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CHOICE_HPP

#include "dca/linalg/linalg.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"
#ifdef DCA_HAVE_CUDA
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_gpu_submatrix.hpp"
#endif  // DCA_HAVE_CUDA

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

namespace {
template <linalg::DeviceType device, class Parameters, bool use_submatrix, typename Real>
struct CtintWalkerChoicheSelector;

template <class Parameters, typename Real>
struct CtintWalkerChoicheSelector<linalg::CPU, Parameters, true, Real> {
  using type = CtintWalkerSubmatrix<linalg::CPU, Parameters, Real>;
};
template <class Parameters, typename Real>
struct CtintWalkerChoicheSelector<linalg::CPU, Parameters, false, Real> {
  using type = CtintWalker<linalg::CPU, Parameters, Real>;
};

template <class Parameters, typename Real>
struct CtintWalkerChoicheSelector<linalg::GPU, Parameters, true, Real> {
  using type = CtintWalkerSubmatrix<linalg::GPU, Parameters, Real>;
};
template <class Parameters, typename Real>
struct CtintWalkerChoicheSelector<linalg::GPU, Parameters, false, Real> {
  // There is only a submatrix implementation of the GPU walker.
};
}  // namespace

template <linalg::DeviceType device_t, class Parameters, bool use_submatrix>
using CtintWalkerChoice =
    typename CtintWalkerChoicheSelector<device_t, Parameters, use_submatrix>::type;

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CHOICE_HPP
