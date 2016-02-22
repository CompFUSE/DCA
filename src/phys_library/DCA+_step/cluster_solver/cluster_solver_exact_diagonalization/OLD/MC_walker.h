//-*-C++-*-

#ifndef MC_WALKER_H
#define MC_WALKER_H

namespace QMC 
{
  /*! 
   * \class   MC_walker
   * \ingroup MONTE-CARLO-INTEGRATOR
   * \brief   empty template for a Monte Carlo walker
   * \author  Peter Staar
   * \version 1.0
   */
  template<MC_integration_method_type MC_integration_method_t, LIN_ALG::device_type device_t, class parameters_type, class base_cluster_type>
  class MC_walker
  {
    typedef typename parameters_type::random_number_generator                                 rng_type;
    typedef MultiOrbitalMultiSiteStructure<parameters_type, base_cluster_type> MOMS_type;

  public:
    
    MC_walker(parameters_type& parameters_ref,
	      MOMS_type&       MOMS_ref,
	      rng_type&        rng_ref);

    ~MC_walker();

    void  initialize();

    bool& is_thermalized();

    void  do_sweep();
    void  do_step();

    int   get_sign();

    template<class stream_type>
    void to_JSON(stream_type& ss);
   };

  /*! 
   * \class   MC_walker
   * \ingroup MONTE-CARLO-INTEGRATOR
   * \brief   empty template for a Monte Carlo BIT-walker (BIT=build-in-test)
   * \author  Peter Staar
   * \version 1.0
   */
  template<MC_integration_method_type MC_integration_method_t, class parameters_type, class base_cluster_type>
  class MC_walker_BIT
  {
    typedef MultiOrbitalMultiSiteStructure<parameters_type, base_cluster_type> MOMS_type;

  public:
    
    MC_walker_BIT(parameters_type& parameters_ref,
		  MOMS_type&       MOMS_ref);

    ~MC_walker_BIT();
  };

}

#endif
