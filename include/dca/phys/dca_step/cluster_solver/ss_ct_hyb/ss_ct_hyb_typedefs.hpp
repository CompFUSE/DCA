// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
//
// This class defines common types for the Single-Site Hybridization Monte Carlo Integrator.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_TYPEDEFS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_TYPEDEFS_HPP

#include "dca/linalg/matrix.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/structures/ss_ct_hyb_configuration.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

template <class parameters_type, class MOMS_type>
class SsCtHybTypedefs {
public:
  // Types that define the profiling.
  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_type;

  // Types that define the scalar type and matrix type.
  typedef double scalartype;
  // typedef resizeable_square_matrix<scalartype> vertex_vertex_matrix_type;
  typedef dca::linalg::Matrix<scalartype, dca::linalg::CPU> vertex_vertex_matrix_type;

  // Types that define the vertex and configuration type.
  typedef SS_CT_HYB_configuration configuration_type;
  typedef typename configuration_type::orbital_configuration_type orbital_configuration_type;
};

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_SS_CT_HYB_TYPEDEFS_HPP
