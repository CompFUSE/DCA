//-*-C++-*-

#ifndef SWAP_SEGMENT_TOOLS_H
#define SWAP_SEGMENT_TOOLS_H

namespace QMC {
  
/*!
 *  \ingroup HYBRIDIZATION
 *
 *  \author  Bart Ydens
 *  \brief   This class implements the swap of segemnts between two hybridization lines.
 *
 */
  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  class swap_segment_tools
  {
#include "type_definitions.h"
   
    typedef resizeable_square_matrix     <double> vertex_vertex_matrix_type;

    typedef typename configuration_t::orbital_configuration_t orbital_configuration_t;
    
    typedef swap_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t> THIS_TYPE;
    
  public:
    
    swap_segment_tools(configuration_t& configuration,
		       parameters_t&    parameters,
		       MOMS_t&          MOMS,
		       concurrency_t&   concurrency);

    ~swap_segment_tools();

    template<typename function_type_0, typename function_type_1, typename function_type_2>
    double swap_orbitals(int i,
		       int j, 
		       function_type_0& mu,
		       double& sign,
		       function_type_1& M,
		       function_type_2& F);

    template<typename function_type_1, typename function_type_2>
    double construct_inverse(function_type_1& M, 
			     double BETA,
			     function_type_2& F, 
			     int flavor_1, 
			     int flavor_2);

    template<typename matrix_type_1, typename function_type_2>
    void construct_matrix(matrix_type_1& M, 
			  double BETA,  
			  function_type_2& F, 
			  int flavor_1,
			  int flavor_2);

  private:
    
    configuration_t& configuration;
        
    parameters_t&    parameters;
    MOMS_t&          MOMS;
    concurrency_t&   concurrency;

    int    FLAVORS;
    double BETA;
  };

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  swap_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::swap_segment_tools(configuration_t& configuration_ref,
											       parameters_t&    parameters_ref,
											       MOMS_t&          MOMS_ref,
											       concurrency_t&   concurrency_ref):
    configuration(configuration_ref),
    
    parameters(parameters_ref),
    MOMS(MOMS_ref),
    concurrency(concurrency_ref)
  {
    FLAVORS = b::dmn_size()*s::dmn_size(); 
    BETA   = parameters.get_beta();
  }

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  swap_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::~swap_segment_tools()
  {}

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename function_type_0, typename function_type_1, typename function_type_2>
  double swap_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::swap_orbitals(int i,
											       int j, 
											       function_type_0& mu,
											       double& sign,
											       function_type_1& M,
											       function_type_2& F)
  {
    orbital_configuration_t& vertices      = configuration.get_vertices(i);
    orbital_configuration_t& swap_vertices = configuration.get_vertices(j);

    vertex_vertex_matrix_type M_new_this, M_new_swap;
  
    double det_old_this = construct_inverse(M_new_this, BETA,  F, i, i); // before swap
    double det_old_swap = construct_inverse(M_new_swap, BETA,  F, j, j); 
    double det_new_this = construct_inverse(M_new_this, BETA,  F, i, j); // after swap
    double det_new_swap = construct_inverse(M_new_swap, BETA,  F, j, i);
    double det_rat = (det_new_this/det_old_this)*(det_new_swap/det_old_swap);
  
    Hybridization_vertex full_segment(0,BETA);

    double length_this = HYBRIDIZATION_TOOLS::compute_overlap(full_segment, vertices     , configuration.get_full_line(i), BETA);
    double length_swap = HYBRIDIZATION_TOOLS::compute_overlap(full_segment, swap_vertices, configuration.get_full_line(j), BETA);
    double e_site      = length_this*(mu(i)-mu(j))+length_swap*(mu(j)-mu(i));  // site energy contribution

    double overlap_u = 0;
    for (int k=0; k<FLAVORS; k++) {
      if (k!=i && k!=j) {
	double overlap=0;
	if(configuration.get_full_line(i)){
	  overlap += HYBRIDIZATION_TOOLS::compute_overlap(full_segment, configuration.get_vertices(k), configuration.get_full_line(k), BETA);
	}
	else
	  for (typename orbital_configuration_t::iterator it=vertices.begin(); it!=vertices.end(); it++) {
	    overlap += HYBRIDIZATION_TOOLS::compute_overlap(*it, configuration.get_vertices(k), configuration.get_full_line(k), BETA);
	  }
	overlap_u += (MOMS.H_interactions(j,k,0)-MOMS.H_interactions(i,k,0))*overlap;//(u(j,k)-u(i,k))*overlap;
	overlap = 0;
	if(configuration.get_full_line(j)){
	  overlap += HYBRIDIZATION_TOOLS::compute_overlap(full_segment, configuration.get_vertices(k), configuration.get_full_line(k), BETA);
	}
	else	  
	  for (typename orbital_configuration_t::iterator it=swap_vertices.begin(); it!=swap_vertices.end(); it++) {
	    overlap += HYBRIDIZATION_TOOLS::compute_overlap(*it, configuration.get_vertices(k), configuration.get_full_line(k), BETA);
	  }
	overlap_u += (MOMS.H_interactions(i,k,0)-MOMS.H_interactions(j,k,0))*overlap;//(u(i,k)-u(j,k))*overlap;
      }
    }

    double log_prob = log(fabs(det_rat)) + (-overlap_u-e_site);
    if (log(concurrency.get_random_number()) < log_prob) {   // length of segments, overlap and phonon part are not changed

      construct_inverse(M(i), BETA,  F, i, j);
      construct_inverse(M(j), BETA,  F, j, i);

      swap(vertices, swap_vertices);
      int dummy1 = configuration.get_full_line(i);
      configuration.get_full_line(i) = configuration.get_full_line(j);
      configuration.get_full_line(j) = dummy1;

      return exp(-overlap_u-e_site);
    }
    return -1;
  }
  
  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename function_type_1, typename function_type_2>
  double swap_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::construct_inverse(function_type_1& M, 
												     double BETA,
												     function_type_2& F, 
												     int flavor_1, 
												     int flavor_2) 
  {
    double det = 1;
    orbital_configuration_t segments = configuration.get_vertices(flavor_2);
   
    if(segments.size() > 0){

      construct_matrix(M, BETA, F, flavor_1, flavor_2);
     
      invert_plan<double> inv_pln(segments.size(), M.get_global_size());
      memcpy(inv_pln.Matrix, &M(0,0), sizeof(double)*M.get_global_size()*M.get_global_size());
      inv_pln.execute_plan();
      memcpy( &M(0,0), inv_pln.inverted_matrix,sizeof(double)*M.get_global_size()*M.get_global_size());

      for(int i = 0; i< M.get_current_size(); i++){
	det *= inv_pln.Matrix[i+M.get_global_size()*i];
      }
    }
    else
      M.resize(0);

    return det;
  }       
  
  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename matrix_type_1, typename function_type_2>
  void swap_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::construct_matrix(matrix_type_1& M, 
												  double BETA,  
												  function_type_2& F, 
												  int flavor_1,
												  int flavor_2) 
  {
    orbital_configuration_t segments = configuration.get_vertices(flavor_2);

    int N = segments.size();

    matrix_type_1 M_new(N,N);

    int row=-1;
    int column=-1;

    static int* coor = new int[2];
    static nu nu_obj;

    nu_obj.linind_2_subind(flavor_1, coor);

    for (typename orbital_configuration_t::iterator it1=segments.begin(); it1!=segments.end(); it1++) {
      row++;
      for (typename orbital_configuration_t::iterator it2=segments.begin(); it2!=segments.end(); it2++) {
	column++;
	
	double argument = it1->t_end()-it2->t_start();
	double sign = 1;
	if (argument<0) {
	  argument += BETA;	  
	  sign = -1;
	}
	M_new(row,column) = HYBRIDIZATION_TOOLS::interpolate_F(coor, argument, F)*sign;
      }
      column = -1;
    }
    M.copy_from(M_new);
  }

}

#endif
