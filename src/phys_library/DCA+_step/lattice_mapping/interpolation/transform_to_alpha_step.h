//-*-C++-*-

#ifndef DCA_TRANSFORM_TO_ALPHA_STEP_H
#define DCA_TRANSFORM_TO_ALPHA_STEP_H

namespace DCA 
{
  /*!
   *  \author Peter Staar
   *  \brief  This class computes the analystical transform, used in the intrpolation of the self-energy.
   *
   *  \f{eqnarray*}{
   *     T     [\Sigma] &=& [ \Sigma - \alpha I ]^{-1}
   *     T^{-1}[\Sigma] &=&  T[\Sigma]^{-1} + \alpha I
   *  \f}
   */
  class transform_to_alpha
  {
#include "type_definitions.h"

  public:

    template<typename scalar_type, typename k_dmn_t>
    void forward(scalar_type alpha, 
		 function<std::complex<scalar_type>, k_dmn_t>&     f_k,
		 function<std::complex<scalar_type>, k_dmn_t>& alpha_k);

    template<typename scalar_type, typename k_dmn_t>
    void backward(scalar_type alpha, 
		  function<std::complex<scalar_type>, k_dmn_t>&     f_k,
		  function<std::complex<scalar_type>, k_dmn_t>& alpha_k);

    template<typename scalar_type, typename k_dmn_t>
    static void forward(scalar_type alpha, 
			function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& f_k_w,
			function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& alpha_k_w);

    template<typename scalar_type, typename k_dmn_t>
    static void backward(scalar_type alpha, 
			 function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& f_k_w,
			 function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& alpha_k_w);

    template<typename scalar_type, typename k_dmn_t, typename w_dmn_t>
    static void forward(scalar_type alpha, 
			function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn_t, w_dmn_t> >& f_k_w,
			function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn_t, w_dmn_t> >& alpha_k_w);

    template<typename scalar_type, typename k_dmn_t, typename w_dmn_t>
    static void backward(scalar_type alpha, 
			 function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn_t, w_dmn_t> >& f_k_w,
			 function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn_t, w_dmn_t> >& alpha_k_w);
  };

  template<typename scalar_type, typename k_dmn_t>
  void transform_to_alpha::forward(scalar_type alpha, 
				   function<std::complex<scalar_type>, k_dmn_t>&     f_k,
				   function<std::complex<scalar_type>, k_dmn_t>& alpha_k)
  {
    std::complex<scalar_type> I(0.,alpha);
  
    for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++)
      alpha_k(k_ind) = 1./(f_k(k_ind)-I);
  }

  template<typename scalar_type, typename k_dmn_t>
  void transform_to_alpha::backward(scalar_type alpha, 
				    function<std::complex<scalar_type>, k_dmn_t>&     f_k,
				    function<std::complex<scalar_type>, k_dmn_t>& alpha_k)
  {
    std::complex<scalar_type> I(0.,alpha);

    for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++)
      f_k(k_ind) = 1./alpha_k(k_ind)+I;
  }

  template<typename scalar_type, typename k_dmn_t>
  void transform_to_alpha::forward(scalar_type alpha, 
				   function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >&     f_k_w,
				   function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& alpha_k_w)
  {
    std::complex<scalar_type> I(0.,alpha);

    int N = nu::dmn_size();

    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU>       f_matrix ("f_matrix" , std::pair<int,int>(N,N));
    LIN_ALG::GEINV<LIN_ALG::CPU>::plan<std::complex<scalar_type> > geinv_obj(f_matrix);

    for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++)
      {
	for(int j=0; j<N; ++j)
	  for(int i=0; i<N; ++i)
	    f_matrix(i,j) = f_k_w(i,j,k_ind);
	
	for(int i=0; i<N; i++)
	  f_matrix(i,i) -= I;

	//LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(f_matrix);
	geinv_obj.execute(f_matrix);

	for(int j=0; j<N; ++j)
	  for(int i=0; i<N; ++i)
	    alpha_k_w(i,j,k_ind) = f_matrix(i,j);	
      }
  }

  template<typename scalar_type, typename k_dmn_t>
  void transform_to_alpha::backward(scalar_type alpha, 
				    function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >&     f_k_w,
				    function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& alpha_k_w)
  {
    std::complex<scalar_type> I(0.,alpha);

    int N = nu::dmn_size();

    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU>       f_matrix ("f_matrix" , std::pair<int,int>(N,N));
    LIN_ALG::GEINV<LIN_ALG::CPU>::plan<std::complex<scalar_type> > geinv_obj(f_matrix);

    for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++)
      {
	for(int j=0; j<N; ++j)
	  for(int i=0; i<N; ++i)
	    f_matrix(i,j) = alpha_k_w(i,j,k_ind);

	//LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(f_matrix);
	geinv_obj.execute(f_matrix);
	
	for(int i=0; i<N; i++)
	  f_matrix(i,i) += I;

	for(int j=0; j<N; ++j)
	  for(int i=0; i<N; ++i)
	    f_k_w(i,j,k_ind) = f_matrix(i,j);	
      }
  }

  template<typename scalar_type, typename k_dmn_t, typename w_dmn_t>
  void transform_to_alpha::forward(scalar_type alpha, 
				   function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn_t, w_dmn_t> >&     f_k_w,
				   function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn_t, w_dmn_t> >& alpha_k_w)
  {
    int N = nu::dmn_size();

    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU>       f_matrix("f_matrix" , std::pair<int,int>(N,N));
    LIN_ALG::GEINV<LIN_ALG::CPU>::plan<std::complex<scalar_type> > geinv_obj(f_matrix);

    for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++)
      {
	scalar_type factor = w_dmn_t::get_elements()[w_ind]>0 ? 1 : -1;

	std::complex<scalar_type> I(0., factor*alpha);

	for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++)
	  {
	    for(int j=0; j<N; ++j)
	      for(int i=0; i<N; ++i)
		f_matrix(i,j) = f_k_w(i,j,k_ind,w_ind);
	    
	    for(int i=0; i<N; i++)
	      f_matrix(i,i) -= I;
	    
	    //LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(f_matrix);
	    geinv_obj.execute(f_matrix);
	    
	    for(int j=0; j<N; ++j)
	      for(int i=0; i<N; ++i)
		alpha_k_w(i,j,k_ind,w_ind) = f_matrix(i,j);	
	  }
      }
  }

  template<typename scalar_type, typename k_dmn_t, typename w_dmn_t>
  void transform_to_alpha::backward(scalar_type alpha, 
				    function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn_t, w_dmn_t> >&     f_k_w,
				    function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn_t, w_dmn_t> >& alpha_k_w)
  {
    int N = nu::dmn_size();

    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU>       f_matrix("f_matrix" , std::pair<int,int>(N,N));
    LIN_ALG::GEINV<LIN_ALG::CPU>::plan<std::complex<scalar_type> > geinv_obj(f_matrix);

    for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++)
      {
	scalar_type factor = w_dmn_t::get_elements()[w_ind]>0 ? 1 : -1;

	std::complex<scalar_type> I(0., factor*alpha);

	for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); k_ind++)
	  {
	    for(int j=0; j<N; ++j)
	      for(int i=0; i<N; ++i)
		f_matrix(i,j) = alpha_k_w(i,j,k_ind,w_ind);
	    
	    //LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(f_matrix);
	    geinv_obj.execute(f_matrix);

	    for(int i=0; i<N; i++)
	      f_matrix(i,i) += I;
	    
	    for(int j=0; j<N; ++j)
	      for(int i=0; i<N; ++i)
		f_k_w(i,j,k_ind,w_ind) = f_matrix(i,j);	
	  }
      }
  }

}

#endif
