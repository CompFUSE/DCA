//-*-C++-*-

#ifndef MAKE_G4_0_MATRIX_H
#define MAKE_G4_0_MATRIX_H

namespace dca {

  /*! \file MAKE_G4_0_MATRIX.h
   *
   *  \author Peter Staar
   */
  template<class parameter_type, class MOMS_type>
  class make_G4_0_matrix 
  {
#include "type_definitions.h"
    
  public:
 
    make_G4_0_matrix(parameter_type& parameters_ref, MOMS_type& MOMS_ref);
    ~make_G4_0_matrix();

    template<typename analysis_t>
    void execute(analysis_t& analysis_ref);

    template<typename scalartype_1, typename scalartype_2, typename k_dmn_t>
    static void execute(parameter_type&                                                                                     parameters,  
			FUNC_LIB::function<scalartype_1, dmn_4<nu,nu,k_dmn_t,w> >&                                            G,
			FUNC_LIB::function<scalartype_2, dmn_2< dmn_4<b,b,k_dmn_t,w_VERTEX>, dmn_4<b,b,k_dmn_t,w_VERTEX> > >& P0);
    

  private:    
    
    parameter_type&  parameters;
    MOMS_type&       MOMS;

    k_DCA_w_VERTEX k_DCA_w_VERTEX_domain;
  };

  template<class parameter_type, class MOMS_type>
  make_G4_0_matrix<parameter_type, MOMS_type>::make_G4_0_matrix (parameter_type& parameters_in, 
								 MOMS_type& MOMS_in):
    parameters(parameters_in),
    MOMS(MOMS_in)
  {}

  template<class parameter_type, class MOMS_type>
  make_G4_0_matrix<parameter_type, MOMS_type>::~make_G4_0_matrix()
  {}

  template<class parameter_type, class MOMS_type>
  template<typename analysis_t>
  void make_G4_0_matrix<parameter_type, MOMS_type>::execute(analysis_t& analysis_ref)
  {
    analysis_ref.G4_0_b_k_w__b_k_w = 0;

    int W        = parameters.get_number_of_positive_frequencies();
    int W_vertex = w_VERTEX::dmn_size()/2;

    int q    = parameters.get_q_channel();
    int w_nu = parameters.get_w_channel();

    int* coor = new int[2];

    for(int i=0; i<k_DCA_w_VERTEX_domain.get_size(); i++) 
      {
	k_DCA_w_VERTEX_domain.linind_2_subind(i, coor);

	int k        = coor[0];
	int w_vertex = coor[1];
	int w        = (coor[1] - W_vertex) + W;
	    
	int k_plus_q  = DCA_k_cluster_type::add     (k, q);     
	int q_minus_k = DCA_k_cluster_type::subtract(k, q); 
	//int min_k = DCA_k_cluster_type::subtract(k, 0); 
	    
	for(int n1=0; n1<b::dmn_size(); n1++){
	  for(int n2=0; n2<b::dmn_size(); n2++){
	    for(int m1=0; m1<b::dmn_size(); m1++){
	      for(int m2=0; m2<b::dmn_size(); m2++){

		switch(parameters.get_vertex_measurement_type()) 
		  {
		  case PARTICLE_HOLE_TRANSVERSE:
		  
		    {
		      double wn           = w::get_elements()[w];
		      double w_nu_plus_wn = w::get_elements()[w+w_nu];
		      double beta         = parameters.get_beta();
		      if(std::fabs( (w_nu*2*M_PI/beta+wn) - w_nu_plus_wn)>1.e-6){
            std::cout << (w_nu*2*M_PI/beta+wn) << "\t" <<  w_nu_plus_wn << std::endl;
			throw std::logic_error(__FUNCTION__);
		      }

		      analysis_ref.G4_0_b_k_w__b_k_w(n1,n2,k,w_vertex,
						     m1,m2,k,w_vertex) = -MOMS.G_k_w(n1, e_UP, m2, e_UP, k, w)*MOMS.G_k_w(n2, e_UP, m1, e_UP, k_plus_q, w+w_nu);
		    }
		    break;

		  case PARTICLE_HOLE_MAGNETIC:
		  
		    analysis_ref.G4_0_b_k_w__b_k_w(n1,n2,k,w_vertex,
						   m1,m2,k,w_vertex) = -MOMS.G_k_w(n1, e_UP, m2, e_UP, k, w)*MOMS.G_k_w(n2, e_UP, m1, e_UP, k_plus_q, w+w_nu);
		    break;
		      
		  case PARTICLE_HOLE_CHARGE:
	
		    analysis_ref.G4_0_b_k_w__b_k_w(n1,n2,k,w_vertex,
						   m1,m2,k,w_vertex) = -2.*MOMS.G_k_w(n1, e_UP, m1, e_UP, k, w)*MOMS.G_k_w(n2, e_UP, m2, e_UP, k_plus_q, w+w_nu);
		    break;
		
		  case PARTICLE_PARTICLE_UP_DOWN:
		    assert(std::fabs(w::get_elements()[w]+w::get_elements()[2*W-1-w]) < 1.e-6);

		    analysis_ref.G4_0_b_k_w__b_k_w(n1,n2,k,w_vertex,
						   m1,m2,k,w_vertex) = MOMS.G_k_w(n1, e_UP, m1, e_UP, k, w)*MOMS.G_k_w(n2, e_UP, m2, e_UP, q_minus_k, w_nu+(2*W-1-w));
// 		    analysis_ref.G4_0_b_k_w__b_k_w(n1,n2,k,w_vertex,
// 						   m1,m2,k,w_vertex) = MOMS.G_k_w(n1, e_UP, m1, e_UP, min_k, 2*W-1-w)*MOMS.G_k_w(n2, e_UP, m2, e_UP, k_plus_q, w);
		    break;
		      
		  default:
		    throw std::logic_error(__FUNCTION__);
		  }
	      }
	    }
	  }
	}
      }
  
    delete [] coor;
  }

  template<class parameter_type, class MOMS_type>
  template<typename scalartype_1, typename scalartype_2, typename k_dmn_t>
  void make_G4_0_matrix<parameter_type, MOMS_type>::execute(parameter_type&                                                                                    parameters,  
							    FUNC_LIB::function<scalartype_1, dmn_4<nu,nu,k_dmn_t,w> >&                                           G,
							    FUNC_LIB::function<scalartype_2, dmn_2< dmn_4<b,b,k_dmn_t,w_VERTEX>,dmn_4<b,b,k_dmn_t,w_VERTEX> > >& P0)
  {
    P0 = 0.;

    dmn_2<k_dmn_t, w_VERTEX> k_w_dmn;

    int W        = parameters.get_sp_fermionic_frequencies();;//parameters.get_number_of_positive_frequencies();
    int W_vertex = w_VERTEX::dmn_size()/2;
    int q        = parameters.get_q_channel();

    int w_nu = parameters.get_w_channel();

    int coor[2];

    for(int i=0; i<k_w_dmn.get_size(); i++) 
      {
	k_w_dmn.linind_2_subind(i, coor);

	int k        = coor[0];
	int w_vertex = coor[1];
	int w        = (coor[1] - W_vertex) + W;
	    
	int k_plus_q  = k_dmn_t::parameter_type::add     (k, q);     
	int q_minus_k = k_dmn_t::parameter_type::subtract(k, q); 
	    
	for(int n1=0; n1<b::dmn_size(); n1++){
	  for(int n2=0; n2<b::dmn_size(); n2++){
	    for(int m1=0; m1<b::dmn_size(); m1++){
	      for(int m2=0; m2<b::dmn_size(); m2++){

		switch(parameters.get_vertex_measurement_type()) 
		  {
		  case PARTICLE_HOLE_TRANSVERSE:
		    {
		      P0(n1,n2,k,w_vertex,
			 m1,m2,k,w_vertex) = -G(n1, e_UP, m2, e_UP, k, w)*G(n2, e_UP, m1, e_UP, k_plus_q, w+w_nu);
		      break;
		    }

		  case PARTICLE_HOLE_MAGNETIC:
		    {		  
		      P0(n1,n2,k,w_vertex,
			 m1,m2,k,w_vertex) = -G(n1, e_UP, m2, e_UP, k, w)*G(n2, e_UP, m1, e_UP, k_plus_q, w+w_nu);
		      break;
		    }
		  case PARTICLE_HOLE_CHARGE:
		    {
		      P0(n1,n2,k,w_vertex,
			 m1,m2,k,w_vertex) = -2.*G(n1, e_UP, m1, e_UP, k, w)*G(n2, e_UP, m2, e_UP, k_plus_q, w+w_nu);
		      break;
		    }

		  case PARTICLE_PARTICLE_UP_DOWN:
		    {	
		      double wn          = w::get_elements()[w];
		      double w_nu_min_wn = w::get_elements()[w_nu+(2*W-1-w)];
		      double beta        = parameters.get_beta();

		      if(std::fabs( (w_nu*M_PI/beta-wn) - w_nu_min_wn)>1.e-6)
            throw std::logic_error(__FUNCTION__);


		      P0(n1,n2,k,w_vertex,
			 m1,m2,k,w_vertex) = G(n1, e_UP, m1, e_UP, k, w)*G(n2, e_UP, m2, e_UP, q_minus_k, w_nu+(2*W-1-w));
		      break;
		    }

		  default:
		    throw std::logic_error(__FUNCTION__);
		  }
	      }
	    }
	  }
	}
      }
  }
}

#endif
