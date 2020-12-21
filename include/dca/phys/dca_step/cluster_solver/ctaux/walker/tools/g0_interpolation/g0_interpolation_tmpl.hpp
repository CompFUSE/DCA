// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.eth.ch)
//
// This class provides a template for the interpolation of
// \f$G^{0}\f$ towards the \f$G^{0}\f$-matrix.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_TMPL_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_TMPL_HPP

#include "dca/linalg/device_type.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <dca::linalg::DeviceType device_t, typename Parameters, typename Real>
class G0Interpolation;

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G0_INTERPOLATION_G0_INTERPOLATION_TMPL_HPP
