//-*-C++-*-

#ifndef DCA_QMCI_ACCUMULATOR_H
#define DCA_QMCI_ACCUMULATOR_H
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_template.h"

namespace DCA {
namespace QMCI {
/*!
 * \class MC_accumulator
 * \brief empty template for a Monte Carlo accumulator
 * \author Peter Staar
 * \version 1.0
 */
template <QMCI_NAMES NAME, LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
class MC_accumulator {
public:
  MC_accumulator(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  ~MC_accumulator();

  template <typename dca_info_struct_t>
  void finalize(dca_info_struct_t& dca_info_struct);

  void initialize();

  void accumulate();

  template <class stream_type>
  void to_JSON(stream_type& ss);
};

/*!
 *  \ingroup CT-AUX
 *
 *  \brief   empty template for the single-particle measurements.
 *  \author  Peter Staar
 *  \version 1.0
 */
template <QMCI_NAMES QMCI_NAME, QMCI_SP_MEASUREMENT_NAMES QMCI_SP_MEASUREMENT_NAME,
          class parameters_type, class MOMS_type>
class MC_single_particle_accumulator {};
}
}

#endif
