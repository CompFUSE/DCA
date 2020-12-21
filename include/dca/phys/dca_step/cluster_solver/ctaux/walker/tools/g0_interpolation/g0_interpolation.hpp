// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.eth.ch)
//
// This file includes all the specializations of the G0Interpolation class.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_HPP

#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_cpu.hpp"
#ifdef DCA_HAVE_CUDA
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_gpu.hpp"
#endif

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_HPP
