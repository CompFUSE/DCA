//-*-C++-*-

#ifndef COMPUTE_LATTICE_GREENS_FUNCTION_H
#define COMPUTE_LATTICE_GREENS_FUNCTION_H

namespace DCA
{

  namespace SERIES_EXPANSION
  {
    template<class parameters_type, class MOMS_type, class k_dmn_t, class w_dmn_t>
    class compute_lattice_Greens_function
    {
#include "type_definitions.h"

      typedef typename parameters_type::concurrency_type concurrency_type;

    public:

      typedef FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu, k_HOST, w> > function_type;

    public:

      compute_lattice_Greens_function(parameters_type& parameters_ref,
				      MOMS_type&       MOMS_ref);

      ~compute_lattice_Greens_function();

      function_type& get_G_k_w()   { return G_k_w;  }
      function_type& get_G0_k_w()  { return G0_k_w; }

      void execute();

    private:

      parameters_type&  parameters;
      concurrency_type& concurrency;
      MOMS_type&        MOMS;

      FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> > G_k_w;
      FUNC_LIB::function<std::complex<double>, dmn_4<nu,nu,k_dmn_t,w_dmn_t> > G0_k_w;
    };

    template<class parameters_type, class MOMS_type, class k_dmn_t, class w_dmn_t>
    compute_lattice_Greens_function<parameters_type, MOMS_type, k_dmn_t, w_dmn_t>::compute_lattice_Greens_function(parameters_type& parameters_ref,
														   MOMS_type&       MOMS_ref):
      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),
      MOMS(MOMS_ref),

      G_k_w ("G_k_w"),
      G0_k_w("G0_k_w")
    {}
    
    template<class parameters_type, class MOMS_type, class k_dmn_t, class w_dmn_t>
    compute_lattice_Greens_function<parameters_type, MOMS_type, k_dmn_t, w_dmn_t>::~compute_lattice_Greens_function()
    {}
    
    template<class parameters_type, class MOMS_type, class k_dmn_t, class w_dmn_t>
    void compute_lattice_Greens_function<parameters_type, MOMS_type, k_dmn_t, w_dmn_t>::execute()
    {
      LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> I_k  ("I_matrix", nu::dmn_size());
      LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> G_inv("G_inv", nu::dmn_size());

      LIN_ALG::GEINV<LIN_ALG::CPU>::plan<std::complex<double> > geinv_obj(G_inv);

      for(int w_ind=0; w_ind<w::dmn_size(); w_ind++)
	{
	  std::complex<double> i_wm_plus_mu;
	  
	  i_wm_plus_mu.real( parameters.get_chemical_potential() );
	  i_wm_plus_mu.imag( w::get_elements()[w_ind] );
	  
	  for(int i=0; i<nu::dmn_size(); i++)
	    I_k(i,i) = i_wm_plus_mu;

	  for(int k_ind=0; k_ind<k_HOST::dmn_size(); k_ind++){
	    
	    {
	      for(int j=0; j<nu::dmn_size(); j++)
		for(int i=0; i<nu::dmn_size(); i++)
		  G_inv(i,j) = I_k(i,j)-MOMS.H_HOST(i,j,k_ind);
	      
	      //LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(G_inv);
	      geinv_obj.execute(G_inv);

	      for(int j=0; j<nu::dmn_size(); j++)
		for(int i=0; i<nu::dmn_size(); i++)
		  G0_k_w(i,j,k_ind,w_ind) = G_inv(i,j);
	    }

	    {
	      for(int j=0; j<nu::dmn_size(); j++)
		for(int i=0; i<nu::dmn_size(); i++)
		  G_inv(i,j) = I_k(i,j)-MOMS.H_HOST(i,j,k_ind)-MOMS.Sigma_lattice(i,j,k_ind,w_ind);
	      
	      //LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(G_inv);
	      geinv_obj.execute(G_inv);

	      for(int j=0; j<nu::dmn_size(); j++)
		for(int i=0; i<nu::dmn_size(); i++)
		  G_k_w(i,j,k_ind,w_ind) = G_inv(i,j);
	    }
	  }
	}     
    }
    
  }

}

#endif
