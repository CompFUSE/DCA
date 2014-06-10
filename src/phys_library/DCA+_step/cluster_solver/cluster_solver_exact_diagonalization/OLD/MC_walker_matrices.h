//-*-C++-*-

#ifndef MC_WALKER_MATRICES_H
#define MC_WALKER_MATRICES_H

namespace QMC 
{
  /*! 
   * \class   MC_walker
   * \ingroup MONTE-CARLO-INTEGRATOR
   * \brief   empty template for the matrices in a Monte Carlo walker.
   * \author  Peter Staar
   * \version 1.0
   */
  template<MC_integration_method_type MC_integration_method_t, LIN_ALG::device_type device_t, typename parameters_type>
  class MC_walker_matrices
  {
  public:
    
    MC_walker_matrices(parameters_type& parameters);

    ~MC_walker_matrices();
   };
}

#endif
