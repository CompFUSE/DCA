//-*-C++-*-

#ifndef ANTI_SEGMENT_TOOLS_H
#define ANTI_SEGMENT_TOOLS_H

namespace QMC {

/*!
 *  \ingroup HYBRIDIZATION
 *
 *  \author  Bart Ydens
 *  \brief   This class implements the insertion and removal of antisegments \f$(c-c^{\dagger} pair)\f$.
 *
 *
 *  \f{eqnarray}{
 *   W_{acc}(c_k \rightarrow c_{k+1}) = min \left(1, \frac{\beta l_{max}}{k+1} \frac{det F(k+1)}{det F(k)}\frac{W_{loc}(k+1)}{W_{loc}(k)} \right)
 *  \f}
 *
 */
  
  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  class anti_segment_tools
  {
#include "type_definitions.h"
    
    typedef typename configuration_t::orbital_configuration_t orbital_configuration_t;

    typedef anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t> THIS_TYPE;

  public:

    anti_segment_tools(configuration_t& configuration,
		       parameters_t&    parameters,
		       MOMS_t&          MOMS,
		       concurrency_t&   concurrency_ref);

    ~anti_segment_tools();
    
    template<typename function_type_1, typename function_type_2>
    double insert_anti_segment(int this_flavor, 
			     double mu,
			     double& sign,
			     function_type_1& M,
			     function_type_2& F);

    template<typename function_type_1, typename function_type_2>
    double remove_anti_segment(int this_flavor, 
			     double mu,
			     double& sign,
			     function_type_1& M,
			     function_type_2& F);

  private:

    double get_other_length_u(int j, Hybridization_vertex& vertex);

    /*******************
     ***  INSERTION  ***
     *******************/

    template<typename function_type_1, typename function_type_2>
    double insert_anti_segment_into_full_line(int this_flavor, 
					    double mu,
					    double& sign,
					    function_type_1& M,
					    function_type_2& F);

    template<typename function_type_1, typename function_type_2>
    double insert_anti_segment_into_segmented_line(int this_flavor, 
						 double mu,
						 double& sign,
						 function_type_1& M,
						 function_type_2& F);

    /*******************
     ***  REMOVAL    ***
     *******************/

    template<typename function_type_1, typename function_type_2>
    double remove_anti_segment_for_one_segment(int this_flavor, 
					     double mu,
					     double& sign,
					     function_type_1& M,
					     function_type_2& F);

    template<typename function_type_1, typename function_type_2>
    double remove_anti_segment_for_multiple_segments(int this_flavor, 
						  double mu,
						  double& sign,
						  function_type_1& M,
						  function_type_2& F);


    /*******************
     ***  WHATEVER   ***
     *******************/

    template<typename Hybridization_function_type, typename vertex_vertex_matrix_type, typename orbital_configuration_type> 
    static double det_rat_insert_anti(int                        this_flavor,
				      Hybridization_vertex&      anti_segment, 
				      vertex_vertex_matrix_type& M, 
				      orbital_configuration_type&   segments_old, 
				      Hybridization_function_type&  F, 
				      double                     BETA, 
				      double &                   det_rat_sign, 
				      double &                   overlap,
				      std::vector<double>&       Q, 
				      std::vector<double>&       R);

    template<typename Hybridization_function_type, typename vertex_vertex_matrix_type, typename orbital_configuration_type> 
    static void compute_M_insert_anti(int                        this_flavor,
				      Hybridization_vertex& anti_segment, 
				      int s, 
				      int r, 
				      vertex_vertex_matrix_type& M, 
				      orbital_configuration_type& segments_old, 
				      Hybridization_function_type& F, 
				      double BETA, 
				      double det_rat, 
				      std::vector<double>& Q,
				      std::vector<double>& R);

    template<typename Hybridization_function_type, typename vertex_vertex_matrix_type, typename orbital_configuration_type> 
    static double det_rat_remove_anti(int this_flavor,
				      Hybridization_vertex anti_segment, 
				      int r, 
				      int s, 
				      vertex_vertex_matrix_type& M, 
				      orbital_configuration_type&   segments_old, 
				      Hybridization_function_type&  F, 
				      double BETA, 
				      double & det_rat_sign); 

    template<class vertex_vertex_matrix_type> 
    static void compute_M_remove_anti(vertex_vertex_matrix_type& M, int s, int r); 

  private:

    configuration_t& configuration;

    concurrency_t&   concurrency;
    parameters_t&    parameters;
    MOMS_t&          MOMS;

    int    FLAVORS;
    double BETA;
  };

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::anti_segment_tools(configuration_t& configuration_ref,
											       parameters_t&    parameters_ref,
											       MOMS_t&          MOMS_ref,
											       concurrency_t&   concurrency_ref):
    configuration(configuration_ref),

    concurrency(concurrency_ref),
    parameters(parameters_ref),
    MOMS(MOMS_ref)
  {
    FLAVORS = b::dmn_size()*s::dmn_size(); 
    BETA   = parameters.get_beta();
  }

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::~anti_segment_tools()
  {}

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename function_type_1, typename function_type_2>
  double anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::insert_anti_segment(int this_flavor, 
												     double mu,
												     double& sign,
												     function_type_1& M,
												     function_type_2& F)
  {
    if (configuration.get_full_line(this_flavor)==true) 
      return insert_anti_segment_into_full_line(this_flavor, mu, sign, M, F);
    else 
      return insert_anti_segment_into_segmented_line(this_flavor, mu, sign, M, F);
  }

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename function_type_1, typename function_type_2>
  double anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::insert_anti_segment_into_full_line(int this_flavor, 
														    double mu,
														    double& sign,
														    function_type_1& M,
														    function_type_2& F)
  {
    double t = BETA*concurrency.get_random_number();
    orbital_configuration_t& vertices = configuration.get_vertices(this_flavor);
    
    double length = HYBRIDIZATION_TOOLS::compute_length(concurrency.get_random_number(), BETA, 0);
    double t_end = (t+length < BETA ? t+length : t+length-BETA);
    
    Hybridization_vertex segment_insert(t_end, t);
    Hybridization_vertex segment_remove(t,t_end);
    
    double log_prob, overlap, det_rat, det_rat_sign;
    std::vector<double> Fs(vertices.size()), Fe(vertices.size());
    det_rat = HYBRIDIZATION_TOOLS::det_rat_up(this_flavor, segment_insert, M(this_flavor), vertices, F, Fs, Fe, det_rat_sign, overlap);
    
    double otherlength_u = get_other_length_u(this_flavor, segment_remove);
    
    log_prob = log(BETA*BETA*det_rat)-length*mu+otherlength_u;
    
    if (log(concurrency.get_random_number()) < log_prob) 
      {
	HYBRIDIZATION_TOOLS::compute_M_up(segment_insert, 0, M(this_flavor), vertices, F, Fs, Fe, det_rat*overlap);
	sign *= det_rat_sign;
	vertices.push_back(segment_insert);
	configuration.get_full_line(this_flavor) = false;
	return exp(-length*mu+otherlength_u);
      }
    return -1;
  }

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename function_type_1, typename function_type_2>
  double anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::insert_anti_segment_into_segmented_line(int this_flavor, 
															 double mu,
															 double& sign,
															 function_type_1& M,
															 function_type_2& F)
  {
    double t = BETA*concurrency.get_random_number();

    double t_up;   // distance to next segment up (t_start)
    double t_down; // distance to next segment down (t_end)

    orbital_configuration_t& vertices = configuration.get_vertices(this_flavor);

    typename orbital_configuration_t::iterator s_up; // iterator of the segment up
    typename orbital_configuration_t::iterator s_down; // iterator of the segment down

    HYBRIDIZATION_TOOLS::compute_intervals(t, BETA, t_up, t_down, vertices, s_up, s_down);
	
    if (t_down<0) 
      { // t does lie on a segment -> it's possible to insert an anti-segment starting from t

	double length = HYBRIDIZATION_TOOLS::compute_length(concurrency.get_random_number(), -t_down, 0);

	Hybridization_vertex segment_shrink(s_down->t_start(),t);
	  
	double t_start = t + length;
	if (t_start > BETA)
	  t_start-=BETA;
		
	Hybridization_vertex segment_insert(t_start, s_down->t_end());
	Hybridization_vertex anti_segment(t,t_start);
	  
	double otherlength_u = get_other_length_u(this_flavor, anti_segment);
	  
	double log_prob, overlap, det_rat, det_rat_sign;
	std::vector<double> Q(vertices.size()),R(vertices.size());
	det_rat = THIS_TYPE::det_rat_insert_anti(this_flavor, anti_segment, M(this_flavor), vertices, F, BETA, det_rat_sign, overlap, Q, R);
	
	log_prob = log(BETA*(-t_down)/(vertices.size()+1)*det_rat)-length*mu+otherlength_u;
	  
	if (log(concurrency.get_random_number()) < log_prob) 
	  {
	    int s, r; // s is the segment which is shifted, r the segment which is inserted
	    s = 0;
	    for (typename orbital_configuration_t::iterator it=vertices.begin(); it!=s_down; it++)
	      s++;
	    if (anti_segment.t_end() > segment_shrink.t_start())
	      r = s+1;
	    else {
	      r = 0;  
	      s++;
	    }

	    THIS_TYPE::compute_M_insert_anti(this_flavor, anti_segment, s, r, M(this_flavor), vertices, F, BETA, det_rat*overlap, Q, R);

	    s_down->set_t_end(t);
	    if (segment_insert.t_start()>vertices.begin()->t_start()) {
	      s_down++;
	      vertices.insert(s_down, segment_insert);
	    }
	    else {
	      vertices.insert(vertices.begin(), segment_insert);	  
	    }
	    return exp(-length*mu+otherlength_u);
	  }
      }
    return -1;
  }

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename function_type_1, typename function_type_2>
  double anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::remove_anti_segment(int this_flavor, 
												     double mu,
												     double& sign,
												     function_type_1& M,
												     function_type_2& F)
  {
    orbital_configuration_t& vertices = configuration.get_vertices(this_flavor);

    switch(vertices.size())
      {
      case 0:
	return -1;
	break;

      case 1:
	return remove_anti_segment_for_one_segment(this_flavor, mu, sign, M, F);
	break;

      default:
	return remove_anti_segment_for_multiple_segments(this_flavor, mu, sign, M, F);
      };

//     if (vertices.size()>1)
//       remove_anti_segment_for_multiple_segment(this_flavor, mu, sign, M, F);
//     else 
//       remove_anti_segment_for_one_segment(this_flavor, mu, sign, M, F);
  }

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename function_type_1, typename function_type_2>
  double anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::remove_anti_segment_for_one_segment(int this_flavor, 
														       double mu,
														       double& sign,
														       function_type_1& M,
														       function_type_2& F)
  {
    orbital_configuration_t& vertices = configuration.get_vertices(this_flavor);
    assert(vertices.size()==1);
    typename orbital_configuration_t::iterator s_down = vertices.begin();
    
    double det_rat = std::fabs(M(this_flavor)(0,0));
    double length = s_down->t_start()-s_down->t_end(); 
    if (length<0) length += BETA;
    Hybridization_vertex anti_segment(s_down->t_end(),s_down->t_start());
    
    double otherlength_u = get_other_length_u(this_flavor, anti_segment);
    
    double log_prob = log(BETA*BETA/det_rat)-length*mu+otherlength_u;
    
    if (log(concurrency.get_random_number()) < -log_prob) { 
      configuration.get_full_line(this_flavor) = true;//full_line=1;
      vertices.erase(s_down);
      HYBRIDIZATION_TOOLS::compute_M_down(0,M(this_flavor)); 
      return exp(length*mu-otherlength_u);
    }  
    return -1;
  }

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename function_type_1, typename function_type_2>
  double anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::remove_anti_segment_for_multiple_segments(int this_flavor, 
															     double mu,
															     double& sign,
															     function_type_1& M,
															     function_type_2& F)
  {
    typename orbital_configuration_t::iterator s_up; // iterator of the segment up
    typename orbital_configuration_t::iterator s_down; // iterator of the segment down

    orbital_configuration_t& vertices = configuration.get_vertices(this_flavor);

    assert(vertices.size()>1); 

    int r = concurrency.get_random_number()*vertices.size();
    s_up = vertices.begin();
    for (int i=0; i<r; i++) s_up++;

    int s = r-1;
    if (s<0) {
      s=vertices.size()-1;
      s_down=vertices.end();
      s_down--;	  
    }
    else {
      s_down = s_up;
      s_down--;
    }
  
    double length = s_up->t_start() - s_down->t_end();
    if (length < 0) length += BETA;
  
    double t_total = s_up->t_end() - s_down->t_end();
    if (t_total < 0) t_total += BETA;	  

    Hybridization_vertex anti_segment(s_down->t_end(),s_up->t_start());
	
    double otherlength_u = get_other_length_u(this_flavor, anti_segment);
	
    double log_prob, det_rat, det_rat_sign;
	
    det_rat = THIS_TYPE::det_rat_remove_anti(this_flavor, anti_segment, r, s, M(this_flavor), vertices, F, BETA, det_rat_sign);
			  
    log_prob = log(BETA*t_total/vertices.size()/det_rat)-length*mu+otherlength_u;
	
    if (log(concurrency.get_random_number()) < -log_prob) 
      {
        THIS_TYPE::compute_M_remove_anti(M(this_flavor), s, r);
	
	double t_end = s_up->t_end();		
	vertices.erase(s_up);
	
	if (r>0) {
	  s_up=vertices.begin();
	  for (int k=0; k<s; k++)
	    s_up++;
	}
	else {
	  s=vertices.size()-1;
	  s_up = vertices.end();
	  s_up--;
	}
	s_up->set_t_end(t_end);
	return exp(length*mu-otherlength_u);
    }
    return -1;
  }


  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  double anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::get_other_length_u(int j, Hybridization_vertex& vertex)
  {
    double otherlength_u=0;

    for (int i=0; i<FLAVORS; i++) 
      {
	if(i==j) 
	  continue;
	    
	double other_length = HYBRIDIZATION_TOOLS::compute_overlap(vertex, configuration.get_vertices(i), configuration.get_full_line(i), BETA);
	otherlength_u += other_length*MOMS.H_interactions(i, j, 0);
      }

    return otherlength_u;
  }


  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename Hybridization_function_type, typename vertex_vertex_matrix_type, typename orbital_configuration_type> 
  double anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::det_rat_insert_anti(int                        this_flavor,
												       Hybridization_vertex&      anti_segment, 
												       vertex_vertex_matrix_type& M, 
												       orbital_configuration_type&   segments_old, 
												       Hybridization_function_type&  F, 
												       double                     BETA, 
												       double &                   det_rat_sign, 
												       double &                   overlap, 
												       std::vector<double>&       Q_prime,
												       std::vector<double>&       R) 
  {
    static gemv_plan  <double> gemv_pln(1, 1);
    static std::vector<double> Q       (1, 0.);
    Q.resize(M.get_current_size());

    static int* coor = new int[2];
    static nu nu_obj;
    nu_obj.linind_2_subind(this_flavor, coor);
    
    typename orbital_configuration_t::iterator it=segments_old.begin();
    for (int i=0; i<M.get_current_size(); i++) {
      R[i] = HYBRIDIZATION_TOOLS::interpolate_F(coor, anti_segment.t_start()-it->t_start(), F);
      Q[i] = HYBRIDIZATION_TOOLS::interpolate_F(coor, it->t_end()-anti_segment.t_end(), F);
      it++;
    }

    HYBRIDIZATION_TOOLS::compute_Q_prime(Q, M, Q_prime);

    double det_rat = -HYBRIDIZATION_TOOLS::interpolate_F(coor, anti_segment.t_start()-anti_segment.t_end(), F);

    for (int i=0; i<M.get_current_size(); i++) {
	det_rat -= R[i]*Q_prime[i];
    }
  
    {
      overlap = 1;
      // take care of sign changes produced by segments which "wind around"
      // check if anti-segment winds around
      if (anti_segment.t_end()<anti_segment.t_start()) {
	det_rat *= -1;
	overlap = -1;	  
      }
      
      if (det_rat < 0) {
	det_rat_sign = -1;
	det_rat *= -1;
      }
      else {
	det_rat_sign = 1;
      }
    }
    return det_rat;
  }

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename Hybridization_function_type, typename vertex_vertex_matrix_type, typename orbital_configuration_type> 
  void anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::compute_M_insert_anti(int                        this_flavor,
												       Hybridization_vertex& anti_segment, 
												       int s, 
												       int r, 
												       vertex_vertex_matrix_type& M, 
												       orbital_configuration_type& segments_old, 
												       Hybridization_function_type& F, 
												       double BETA, 
												       double S_prime_inv, 
												       std::vector<double>& Q_prime,
												       std::vector<double>& R) 
  {
    static gemv_plan  <double> gemv_pln(1, 1);
    static gemm_plan  <double> gemm_pln;
    static std::vector<double> R_prime (1, 0.);
    if(M.get_current_size()>(int)R_prime.size())
      R_prime.resize(M.get_current_size());

    int i_new, j_new;
    int size=M.get_current_size();

    HYBRIDIZATION_TOOLS::compute_R_prime(R, M, R_prime);
    HYBRIDIZATION_TOOLS::compute_M(Q_prime, R_prime,-1./S_prime_inv, M);

    if(r==0) // need to permute indices of R, L, M
      M.cycle_column_forward();	

    if(r < M.get_current_size() && s < M.get_current_size())
      M.insert_row_and_column(r,s);
    else if(r < M.get_current_size() && s == M.get_current_size())
      M.insert_row_and_add_column(r);
    else if(r == M.get_current_size() && s < M.get_current_size())
      M.add_row_and_insert_column(s);
    else
      M.resize(M.get_current_size()+1);

    M(r,s) = -1./S_prime_inv;

    // row k+1 and column k
    if (r!=0) { // segments remain in the usual order
      for (int i=0; i<size; i++) {
	i_new = (i<r ? i : i+1);
	j_new = (i<s ? i : i+1);
	
	M(i_new,s) = -Q_prime[i]/S_prime_inv;
	M(r,j_new) = -R_prime[i]/S_prime_inv;
      }
    }
    else { // need to permute indices of R, L, M
      for (int i=0; i<size; i++) {
	i_new = (i<r ? i : i+1);
	j_new = (i<s ? i : i+1);
	
	M(i_new,s) = -Q_prime[i]/S_prime_inv;
	M(r,j_new) = -R_prime[HYBRIDIZATION_TOOLS::cycle(i,size)]/S_prime_inv;
      }
    }
  }
  
  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename Hybridization_function_type, typename vertex_vertex_matrix_type, typename orbital_configuration_type>
  double anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::det_rat_remove_anti(int this_flavor,
												       Hybridization_vertex anti_segment, 
												       int r, 
												       int s, 
												       vertex_vertex_matrix_type& M, 
												       orbital_configuration_type&  segments_old, 
												       Hybridization_function_type& F, 
												       double BETA, 
												       double& det_rat_sign) 
  {
    // r is the index of the segment which is removed
    // s is the index of the segment which is shifted

    typename orbital_configuration_type::iterator it=segments_old.begin();
    typename orbital_configuration_type::iterator its(it), itr(it);
    advance(its, s); 
    advance(itr, r);
  
    static int* coor = new int[2];
    static nu nu_obj;
    nu_obj.linind_2_subind(this_flavor, coor);

    double inv_det_rat = -HYBRIDIZATION_TOOLS::interpolate_F(coor, its->t_end()-itr->t_start(), F);

    for (size_t i=0; i<segments_old.size(); i++) {
      if (int(i)!=s) {
	inv_det_rat -= HYBRIDIZATION_TOOLS::interpolate_F(coor, it->t_end()-itr->t_start(), F)*M(r,i)/M(r,s);
      }
      it++;
    }
  
    // take care of sign changes produced by segments which "wind around"
    if (anti_segment.t_end() < anti_segment.t_start()) {
      inv_det_rat *= -1;
    }
	
    if (inv_det_rat < 0) {
      det_rat_sign = -1;
      inv_det_rat *= -1;
    }
    else {
      det_rat_sign = 1;
    }
  
    return 1./inv_det_rat;
  }

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename vertex_vertex_matrix_type> 
  void anti_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::compute_M_remove_anti(vertex_vertex_matrix_type& M, int s, int r) 
  {
    static gemm_plan  <double> gemm_pln;
    static std::vector<double> Q_prime (1, 0.);
    static std::vector<double> R_prime (1, 0.);

    Q_prime.resize(M.get_current_size());
    R_prime.resize(M.get_current_size());

    for(int i=0; i<M.get_current_size(); i++) {
      Q_prime[i] = M(i,s);
      R_prime[i] = M(r,i);
    } 

    HYBRIDIZATION_TOOLS::compute_M(Q_prime, R_prime,-1./Q_prime[r], M);
    
    {      
      M.remove_row_and_column(r,s);

      if(r==0) // need to permute indices of R, L, M
	M.cycle_column_backward();
    }    
  }


} // QMC-namespace
#endif
