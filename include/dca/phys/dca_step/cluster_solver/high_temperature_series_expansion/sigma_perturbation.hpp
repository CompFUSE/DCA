// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Perturbation expansion of the self-energy up to 4th order.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_SIGMA_PERTURBATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_SIGMA_PERTURBATION_HPP

#include <complex>
#include <iostream>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/parallel/util/get_bounds.hpp"
#include "dca/parallel/thread_manager_sum.hpp"
#include "dca/phys/dca_step/cluster_solver/high_temperature_series_expansion/compute_bubble.hpp"
#include "dca/phys/dca_step/cluster_solver/high_temperature_series_expansion/compute_interaction.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace htseries {
// dca::phys::solver::htseries::

// Empty class template.
// N = order of the correction term.
template <int N, class parameter_type, class k_dmn_t>
class sigma_perturbation {};

// Specialization for computation of 1st order term.
#include "sigma_perturbation_1st_order.inc"

// Specialization for computation of 2nd order term.
#include "sigma_perturbation_2nd_order.inc"

// Specialization for computation of 3rd order term.
#include "sigma_perturbation_3rd_order.inc"

// Specialization for computation of 4th order term.
#include "sigma_perturbation_4th_order.inc"

}  // htseries
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_SIGMA_PERTURBATION_HPP
