//-*-C++-*-

#ifndef MC_INTEGRATOR_H
#define MC_INTEGRATOR_H

namespace DCA
{
namespace QMC 
{
  /*!
   *  \defgroup MONTE-CARLO-INTEGRATOR
   */

  /*! 
   * \class   MC_integrator
   * \ingroup MONTE-CARLO-INTEGRATOR
   * \brief   empty template for a Monte Carlo integrator
   * \author  Peter Staar
   * \version 1.0
   */
  template<MC_integration_method_type MC_integration_method_t, LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
  class MC_integrator
  {

  public:
    
    MC_integrator(parameters_type&   parameters_ref,
		  MOMS_type&         MOMS_ref);
    
    ~MC_integrator();

    void initialize(int dca_iteration); 

    void execute();

    template<typename dca_info_struct_t>
    void finalize(dca_info_struct_t& dca_info_struct); 

  };
}

}

#endif
