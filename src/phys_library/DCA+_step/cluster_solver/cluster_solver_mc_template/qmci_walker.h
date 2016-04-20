//-*-C++-*-

#ifndef DCA_QMCI_WALKER_H
#define DCA_QMCI_WALKER_H

namespace DCA {
namespace QMCI {
/*!
 * \class   MC_walker
 * \ingroup MONTE-CARLO-INTEGRATOR
 * \brief   empty template for a Monte Carlo walker
 * \author  Peter Staar
 * \version 1.0
 */
template <QMCI_NAMES NAME, LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
class MC_walker {
  typedef typename parameters_type::random_number_generator rng_type;

  // typedef MultiOrbitalMultiSiteStructure<parameters_type> MOMS_type;

public:
  MC_walker(parameters_type& parameters_ref, MOMS_type& MOMS_ref, rng_type& rng_ref);

  ~MC_walker();

  void initialize();

  bool& is_thermalized();

  void do_sweep();
  void do_step();

  int get_sign();

  template <IO::FORMAT DATA_FORMAT>
  void read(IO::reader<DATA_FORMAT>& reader);

  template <IO::FORMAT DATA_FORMAT>
  void write(IO::writer<DATA_FORMAT>& reader);
};

/*!
 * \class   MC_walker
 * \ingroup MONTE-CARLO-INTEGRATOR
 * \brief   empty template for a Monte Carlo BIT-walker (BIT=build-in-test)
 * \author  Peter Staar
 * \version 1.0
 */
template <QMCI_NAMES NAME, class parameters_type, class MOMS_type>
class MC_walker_BIT {
  // typedef MultiOrbitalMultiSiteStructure<parameters_type> MOMS_type;

public:
  MC_walker_BIT(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  ~MC_walker_BIT();
};
}
}

#endif
