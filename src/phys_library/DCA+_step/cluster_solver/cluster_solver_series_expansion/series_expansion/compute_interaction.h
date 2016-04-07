//-*-C++-*-

#ifndef COMPUTE_INTERACTION_H
#define COMPUTE_INTERACTION_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace DCA
{
  namespace SERIES_EXPANSION 
  {
    /*!
     * \class  compute_interaction
     *
     * \author Peter Staar
     * 
     * \brief  This class implements the interaction matrix
     */
    class compute_interaction
    {

    public:

      typedef FUNC_LIB::function<double, dmn_2<nu,nu> > function_type;

    public:

      compute_interaction()
      {}
     
      ~compute_interaction()
      {}

      template<class r_dmn_t>
      void execute(FUNC_LIB::function<double, dmn_3<nu,nu,r_dmn_t> >& H_interation)
      {
	for(int nu0=0; nu0<2*b::dmn_size(); ++nu0)
	  for(int nu1=0; nu1<2*b::dmn_size(); ++nu1)
	    U(nu0,nu1) = H_interation(nu0,nu1,0);

// 	for(int s_0=0; s_0<2*b::dmn_size(); ++s_0){
// 	  for(int s_1=0; s_1<2*b::dmn_size(); ++s_1)
// 	    cout << U(s_0, s_1) << "\t";
// 	  cout << "\n";
// 	}
// 	cout << "\n";
      }

      function_type& get_function()
      {
	return U;
      }

    protected:

      FUNC_LIB::function<double, dmn_2<nu,nu> > U;  
    };

  }

}

#endif
