//-*-C++-*-

#ifndef MC_ACCUMULATOR_H
#define MC_ACCUMULATOR_H

namespace QMC {

  /*! 
   * \class MC_accumulator
   * \brief empty template for a Monte Carlo accumulator
   * \author Peter Staar
   * \version 1.0
   */
//   template<MC_integration_method_type MC_integration_method_t, class parameters_type, class MOMS_type, class base_cluster_type>//, MC_accumulator_method_type MC_accumulator_method_t = MATSUBARA>
//   class MC_accumulator
//   {

//   public:
    
//     MC_accumulator(parameters_type&   parameters_ref,
// 		   MOMS_type&         MOMS_ref);

//     ~MC_accumulator();

//     template<typename dca_info_struct_t>
//     void finalize(dca_info_struct_t& dca_info_struct); 

//     void initialize(); 

//     void accumulate();

//     template<class stream_type>
//     void to_JSON(stream_type& ss);
//   };
  template<MC_integration_method_type MC_integration_method_t, LIN_ALG::device_type device_t, class parameters_type, class base_cluster_type>
  class MC_accumulator
  {
    typedef MultiOrbitalMultiSiteStructure<parameters_type, base_cluster_type> MOMS_type;

  public:
    
    MC_accumulator(parameters_type&   parameters_ref,
		   MOMS_type&         MOMS_ref);

    ~MC_accumulator();

    template<typename dca_info_struct_t>
    void finalize(dca_info_struct_t& dca_info_struct); 

    void initialize(); 

    void accumulate();

    template<class stream_type>
    void to_JSON(stream_type& ss);
  };

  /*! 
   *  \ingroup CT-AUX
   *
   *  \brief   empty template for the single-particle measurements.
   *  \author  Peter Staar
   *  \version 1.0
   */
  template<MC_integration_method_type MC_integration_method, 
	   MC_accumulator_method_type MC_accumulator_method, 
	   class parameters_type,
	   class base_cluster_type>
  class MC_single_particle_accumulator
  {};
}

#endif
