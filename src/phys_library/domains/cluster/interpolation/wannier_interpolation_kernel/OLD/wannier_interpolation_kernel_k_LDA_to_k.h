//-*-C++-*-

#ifndef WANNIER_INTERPOLATION_KERNEL_k_LDA_TO_k_H
#define WANNIER_INTERPOLATION_KERNEL_k_LDA_TO_k_H

/*! \class   wannier_interpolation_kernel
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   This class implements a Wannier interpolation technique.
 */
template<typename source_dmn_type, typename target_dmn_type>
class wannier_interpolation_kernel
{
#include "type_definitions.h"

//   typedef typename source_dmn_type::parameter_type source_parameter_type;
//   typedef typename target_dmn_type::parameter_type target_parameter_type;

  typedef source_dmn_type source_parameter_type;
  typedef target_dmn_type target_parameter_type;

  typedef typename source_parameter_type::base_cluster source_base_cluster_type;

  typedef r_cluster<FULL, source_base_cluster_type> source_r_cluster_cluster_type;
  typedef k_cluster<FULL, source_base_cluster_type> source_k_cluster_cluster_type;
  
  typedef dmn_0<source_r_cluster_cluster_type> source_r_dmn_t;

  typedef dmn_0<source_k_cluster_cluster_type> source_k_dmn_t;
  typedef dmn_0<target_parameter_type>         target_k_dmn_t;

public:

  wannier_interpolation_kernel();
  ~wannier_interpolation_kernel();

  void execute(std::complex<double>* input_ptr,
	       std::complex<double>* output_ptr);

  FUNC_LIB::function<std::complex<double>, dmn_2<target_k_dmn_t, source_k_dmn_t> >& get_kernel() {return kernel;}

  void reset_kernel();

private:

  void find_phase_factors(bool REDO_INITIALIZE=false);

private:

  function <std::complex<double>, dmn_2<target_k_dmn_t, source_k_dmn_t> > kernel; 
  gemv_plan<std::complex<double> >                                        gemv_pln;
};

template<typename source_dmn_type, typename target_dmn_type>
wannier_interpolation_kernel<source_dmn_type, target_dmn_type>::wannier_interpolation_kernel():
  kernel("wannier_kernel_function"),
  gemv_pln(target_dmn_type::get_size(), source_dmn_type::get_size())
{
  find_phase_factors();

//   for(int K_ind0=0; K_ind0<target_k_dmn_t::dmn_size(); K_ind0++){
//     for(int K_ind1=0; K_ind1<source_k_dmn_t::dmn_size(); K_ind1++){
//       cout << real(kernel( K_ind0,  K_ind1)) << "\t";
//     }
//     cout << endl;
//   }
//   cout << endl;
}

template<typename source_dmn_type, typename target_dmn_type>
wannier_interpolation_kernel<source_dmn_type, target_dmn_type>::~wannier_interpolation_kernel()
{}

template<typename source_dmn_type, typename target_dmn_type>
void wannier_interpolation_kernel<source_dmn_type, target_dmn_type>::execute(std::complex<double>* input_ptr,
									     std::complex<double>* output_ptr)
{
  gemv_pln.vector_source = &input_ptr[0];
  gemv_pln.vector_target = &output_ptr[0];
  gemv_pln.matrix        = &kernel(0);

  gemv_pln.execute_plan();

  gemv_pln.vector_source = NULL;
  gemv_pln.vector_target = NULL;
  gemv_pln.matrix        = NULL;
}

template<typename source_dmn_type, typename target_dmn_type>
void wannier_interpolation_kernel<source_dmn_type, target_dmn_type>::reset_kernel()
{
  kernel.reset();

  find_phase_factors(true);
}

template<typename source_dmn_type, typename target_dmn_type>
void wannier_interpolation_kernel<source_dmn_type, target_dmn_type>::find_phase_factors(bool REDO_INITIALIZE)
{
  {
    FUNC_LIB::function<std::complex<double>, dmn_2<source_r_dmn_t, source_k_dmn_t> > exp_iKR;
    {
      for(int R_ind=0; R_ind<source_r_dmn_t::dmn_size(); R_ind++){
	
	std::vector<double> R_vec = source_r_dmn_t::get_elements()[R_ind];
	
	for(int K_ind=0; K_ind<source_k_dmn_t::dmn_size(); K_ind++){
	  
	  std::vector<double> K_vec = source_k_dmn_t::get_elements()[K_ind];
	  
	  double dot_product=0;
	  for(size_t d=0; d<R_vec.size(); d++)
	    dot_product += R_vec[d]*K_vec[d];
	    
	  exp_iKR(R_ind,K_ind) = std::complex<double>(cos(dot_product), sin(dot_product))/double(source_r_dmn_t::dmn_size());
	}
      }
    }

    FUNC_LIB::function<std::complex<double>, dmn_2<target_k_dmn_t, source_r_dmn_t> > exp_min_ikR;
    {
      FUNC_LIB::function<std::vector<std::vector<double> >, source_r_dmn_t> centered_source_r;

      for(int R_ind=0; R_ind<source_r_dmn_t::dmn_size(); R_ind++)
	centered_source_r(R_ind) = source_r_cluster_cluster_type::find_equivalent_vectors_with_minimal_distance_to_origin(source_r_dmn_t::get_elements()[R_ind]).second;
	
      for(int K_ind=0; K_ind<target_k_dmn_t::dmn_size(); K_ind++){
	  
	std::vector<double> K_vec = target_k_dmn_t::get_elements()[K_ind];
	  
	for(int R_ind=0; R_ind<source_r_dmn_t::dmn_size(); R_ind++){
	    
	  std::complex<double> phase_factor=0;
	    
	  for(size_t r_ind=0; r_ind<centered_source_r(R_ind).size(); r_ind++){
	      
	    std::vector<double> R_vec = centered_source_r(R_ind)[r_ind];
	      
	    double dot_product=0;
	    for(size_t d=0; d<R_vec.size(); d++)
	      dot_product -= R_vec[d]*K_vec[d];
	      
	    phase_factor += std::complex<double>(cos(dot_product), sin(dot_product));
	  }
	    
	  exp_min_ikR(K_ind, R_ind) = phase_factor/double(centered_source_r(R_ind).size());	  
	}
      }
    }

    for(int K_ind0=0; K_ind0<target_k_dmn_t::dmn_size(); K_ind0++)
      for(int K_ind1=0; K_ind1<source_k_dmn_t::dmn_size(); K_ind1++)
	for(int R_ind=0; R_ind<source_r_dmn_t::dmn_size(); R_ind++)
	  kernel(K_ind0, K_ind1) += exp_min_ikR(K_ind0,R_ind)*exp_iKR(R_ind,K_ind1);

    for(int K_ind0=0; K_ind0<target_k_dmn_t::dmn_size(); K_ind0++)
      for(int K_ind1=0; K_ind1<source_k_dmn_t::dmn_size(); K_ind1++)
	if(abs(kernel(K_ind0, K_ind1))<1.e-10)
	  kernel(K_ind0, K_ind1) = 0.;
  }
}

#endif
