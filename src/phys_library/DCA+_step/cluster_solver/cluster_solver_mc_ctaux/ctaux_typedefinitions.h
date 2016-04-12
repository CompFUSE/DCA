//-*-C++-*-

#ifndef DCA_QMCI_CTAUX_TYPE_DEFINITIONS_H
#define DCA_QMCI_CTAUX_TYPE_DEFINITIONS_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_template/qmci_type_definitions.h"

namespace DCA {
namespace QMCI {
/*!
 *
 * \brief   This class defines common types for the CT-AUX Monte Carlo integrator
 * \author  Peter Staar
 * \version 1.0
 */
template <class parameters_type, class MOMS_type>
class MC_type_definitions<CT_AUX_SOLVER, parameters_type, MOMS_type> {
public:
  /*!
   * \brief types that define the profiling
   */
  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_type;
};
}
}

#endif
