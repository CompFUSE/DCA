//-*-C++-*-

#ifndef FULL_LINE_TOOLS_H
#define FULL_LINE_TOOLS_H

namespace QMC {
  /*!
 *  \ingroup HYBRIDIZATION
 *
 *  \author  Bart Ydens
 *  \brief   This class implements the insertion and removal of a full line. The insertion of a full line is only possible if no segement is present on that line.
 *
 *  \f{eqnarray}{
 *   W_{acc}(empty \rightarrow full) = min \left(1, \frac{W_{loc}(full)}{W_{loc}(empty)} \right)
 *  \f}
 *
 */
  template<typename configuration_type, typename parameters_type, typename MOMS_type, typename rng_type>
  class full_line_tools
  {
#include "type_definitions.h"
   
    typedef typename configuration_type::orbital_configuration_type orbital_configuration_type;
    
    typedef full_line_tools<configuration_type, parameters_type, MOMS_type, rng_type> THIS_TYPE;
    
  public:
    
    full_line_tools(configuration_type& configuration,
		    parameters_type&    parameters,
		    MOMS_type&          MOMS,
		    rng_type&           rng_ref);

    ~full_line_tools();

    void insert_full_line(int j, double mu);
    void remove_full_line(int j, double mu);

  private:

    double get_other_length_u(int j);

  private:
    
    configuration_type& configuration;
        
    parameters_type&    parameters;
    MOMS_type&          MOMS;
    rng_type&           rng;

    int    FLAVORS;
    double BETA;
  };

  template<typename configuration_type, typename parameters_type, typename MOMS_type, typename rng_type>
  full_line_tools<configuration_type, parameters_type, MOMS_type, rng_type>::full_line_tools(configuration_type& configuration_ref,
											     parameters_type&    parameters_ref,
											     MOMS_type&          MOMS_ref,
											     rng_type&           rng_ref):
    configuration(configuration_ref),
    
    parameters(parameters_ref),
    MOMS(MOMS_ref),
    rng(rng_ref)
  {
    FLAVORS = b::dmn_size()*s::dmn_size(); 
    BETA   = parameters.get_beta();
  }

  template<typename configuration_type, typename parameters_type, typename MOMS_type, typename rng_type>
  full_line_tools<configuration_type, parameters_type, MOMS_type, rng_type>::~full_line_tools()
  {}

  template<typename configuration_type, typename parameters_type, typename MOMS_type, typename rng_type>
  void full_line_tools<configuration_type, parameters_type, MOMS_type, rng_type>::insert_full_line(int j, double mu)
  {
//     cout << __FUNCTION__ << endl;

    if(configuration.get_full_line(j)==true) 
      return;
  
    double otherlength_u = get_other_length_u(j);
    
    if (log(rng.get_random_number()) < BETA*mu-otherlength_u)
      configuration.get_full_line(j) = true;
  }

  template<typename configuration_type, typename parameters_type, typename MOMS_type, typename rng_type>
  void full_line_tools<configuration_type, parameters_type, MOMS_type, rng_type>::remove_full_line(int j, double mu)
  {
//     cout << __FUNCTION__ << endl;

    if(configuration.get_full_line(j)==false) 
      return;

    double otherlength_u = get_other_length_u(j);

    if (log(rng.get_random_number()) < -BETA*mu+otherlength_u)
      configuration.get_full_line(j) = false;  
  }

  template<typename configuration_type, typename parameters_type, typename MOMS_type, typename rng_type>
  double full_line_tools<configuration_type, parameters_type, MOMS_type, rng_type>::get_other_length_u(int j)
  {
    double otherlength_u=0;

    for(int i=0; i<FLAVORS; i++) 
      {
	if(i==j) 
	  continue;
      
	double other_length=0;
	
	if(configuration.get_full_line(i) == true)
	  other_length = BETA;
	else
	  {
	    orbital_configuration_type& vertices = configuration.get_vertices(i);

	    for(size_t l=0; l<vertices.size(); l++)
	      other_length += (vertices[l].t_end()-vertices[l].t_start()>0 ? vertices[l].t_end()-vertices[l].t_start() : vertices[l].t_end()-vertices[l].t_start()+BETA);
	  }
	
	otherlength_u += other_length*MOMS.H_interactions(i, j, 0);
      }
    
    return otherlength_u;
  }

}

#endif
