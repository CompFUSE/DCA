// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Selects the appropriate ctint walker implementation.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CHOICE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CHOICE_HPP


#include "dca/linalg/linalg.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_gpu_submatrix.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

namespace {
template <linalg::DeviceType device, class Parameters, bool use_submatrix>
struct CtintWalkerChoicheSelector;

template <class Parameters>
struct CtintWalkerChoicheSelector<linalg::CPU, Parameters, true> {
  using type = CtintWalkerSubmatrix<linalg::CPU, Parameters>;
};
template <class Parameters>
struct CtintWalkerChoicheSelector<linalg::CPU, Parameters, false> {
  using type = CtintWalker<linalg::CPU, Parameters>;
};

template <class Parameters>
struct CtintWalkerChoicheSelector<linalg::GPU, Parameters, true> {
  using type = CtintWalkerSubmatrix<linalg::GPU, Parameters>;
};
template <class Parameters>
struct CtintWalkerChoicheSelector<linalg::GPU, Parameters, false> {
  // There is only a submatrix implementation of th GPU walker.
};
}

template <linalg::DeviceType device_t, class Parameters, bool use_submatrix>
using CtintWalkerChoice =
    typename CtintWalkerChoicheSelector<device_t, Parameters, use_submatrix>::type;

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_CTINT_WALKER_CHOICE_HPP
