//-*-C++-*-

#ifndef SHIFT_SEGMENT_TOOLS_H
#define SHIFT_SEGMENT_TOOLS_H

namespace QMC {
  
/*!
 *  \ingroup HYBRIDIZATION
 *
 *  \author  Bart Ydens
 *  \brief   This class implements the shifting of a \f$c^{\dagger}\f$.
 *
 *
 *  \f{eqnarray}{
 *   W_{acc}(c_k \rightarrow c'_k) = min \left(1, \frac{det F'}{det F}\frac{W'_{loc}}{W_{loc}} \right)
 *  \f}
 *
 */

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  class shift_segment_tools
  {
#include "type_definitions.h"
   
    typedef typename configuration_t::orbital_configuration_t orbital_configuration_t;
    
    typedef shift_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t> THIS_TYPE;
    
  public:
    
    shift_segment_tools(configuration_t& configuration,
			parameters_t&    parameters,
			MOMS_t&          MOMS,
			concurrency_t&   concurrency);

    ~shift_segment_tools();

    template<typename function_type_1, typename function_type_2>
    double shift_segment(int this_flavor, 
		       double mu,
		       double& sign,
		       function_type_1& M,
		       function_type_2& F);

  private:

    template<typename vertex_vertex_matrix_t, typename hybridization_function_t>
    double det_rat_shift(int this_flavor,
			 Hybridization_vertex& new_segment, 
			 int k, 
			 vertex_vertex_matrix_t& M,
			 orbital_configuration_t& segments_old, 
			 hybridization_function_t& F, 
			 std::vector<double>&      R,
			 std::vector<double>&      Q_prime,  
			 double & det_rat_sign, 
			 double & overlap); 

    template<typename vertex_vertex_matrix_t, typename hybridization_function_t>
    void compute_M_shift(int this_flavor,
			 Hybridization_vertex& new_segment, 
			 int k, 
			 vertex_vertex_matrix_t& M,
			 orbital_configuration_t& segments_old, 
			 hybridization_function_t& F,
			 std::vector<double>&      R,
			 std::vector<double>&      Q_prime,  
			 double det_rat);

  private:
    
    configuration_t& configuration;
        
    parameters_t&    parameters;
    MOMS_t&          MOMS;
    concurrency_t&   concurrency;

    int    FLAVORS;
    double BETA;
  };

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  shift_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::shift_segment_tools(configuration_t& configuration_ref,
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
  shift_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::~shift_segment_tools()
  {}

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename function_type_1, typename function_type_2>
  double shift_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::shift_segment(int this_flavor, 
												double mu,
												double& sign,
												function_type_1& M,
												function_type_2& F)
  {
    orbital_configuration_t& vertices = configuration.get_vertices(this_flavor);

    int size = vertices.size();

    if (size < 1) 
      return -1;

    int n = size*concurrency.get_random_number();
  
    typename orbital_configuration_t::iterator s, s_up;
    s=vertices.begin();
    for (int i=0; i<n; i++) s++;
    s_up = s; s_up++;
    if (s_up == vertices.end()) 
      s_up = vertices.begin();

    double interval = s_up->t_start() - s->t_start();
    if (interval <= 0) 
      interval += BETA;
  
    double length = HYBRIDIZATION_TOOLS::compute_length(concurrency.get_random_number(), interval, 0);
    double length_old = s->t_end()-s->t_start();
    if (length_old<0)
      length_old += BETA;
  
    double new_t_end = s->t_start() + length;
    if (new_t_end > BETA)
      new_t_end -= BETA;

    Hybridization_vertex segment_insert(s->t_start(), new_t_end);
    Hybridization_vertex segment_remove=*s;

    double otherlength_u=0;
    for (int i=0; i<FLAVORS; i++) 
      {
	if (i==this_flavor) 
	  continue;

	double other_length = HYBRIDIZATION_TOOLS::compute_overlap(segment_insert, configuration.get_vertices(i), configuration.get_full_line(i), BETA)
	  - HYBRIDIZATION_TOOLS::compute_overlap(segment_remove, configuration.get_vertices(i), configuration.get_full_line(i), BETA);

	otherlength_u += other_length*MOMS.H_interactions(i, this_flavor, 0);
    }

    double det_rat, det_rat_sign, overlap;
    std::vector<double> R(vertices.size()), Q(vertices.size());

    det_rat = det_rat_shift(this_flavor, segment_insert, n, M(this_flavor), vertices, F, R, Q, det_rat_sign, overlap);
	
    if (log(concurrency.get_random_number()) < log(det_rat)+(length-length_old)*mu-otherlength_u) 
      {
	compute_M_shift(this_flavor, segment_insert, n, M(this_flavor), vertices, F, R, Q, det_rat*overlap);	
	sign *= det_rat_sign;
	s->set_t_end(new_t_end);
	return exp((length-length_old)*mu-otherlength_u);
      }
    return -1;
  }

  /*! 
   * \ingroup  HYBRIDIZATION
   *
   * \brief    Calculates the determinant ratio for shifting a vertex end point. (u is k-th unity vector)
   * \f{eqnarray}{
   * det_rat &=& 1 + v*A^{-1}*u\\
   *         &=& 1 + (F_{new} - F_{old}) * A^{-1} *u\\
   *         &=& 1 + F_{new} * A^{-1} *u -  F_{old} * A^{-1} *u
   * \f}
   * \f$ F_{old} \f$ is k-th row of matrix A, and \f$A^{-1} *u\f$ is k_th column of \f$A^{-1}\f$ and thus \f$ F_{old} *A^{-1} *u\f$ is equal to 1. (\f$A A^{-1} = I\f$)
   *
   */
  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename vertex_vertex_matrix_t, typename hybridization_function_t>
  double shift_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::det_rat_shift(int this_flavor,
												  Hybridization_vertex& new_segment, 
												  int k, 
												  vertex_vertex_matrix_t& M,
												  orbital_configuration_t& segments_old, 
												  hybridization_function_t& F, 
												  std::vector<double>&      R, 
												  std::vector<double>&      Q_prime,
												  double & det_rat_sign, 
												  double & overlap) 
  {
    static int inc = 1;
    static int* coor = new int[2];
    static nu nu_obj;
    nu_obj.linind_2_subind(this_flavor, coor);

    BLAS::dcopy_(&M.get_current_size(), &M(0,k), &inc, &Q_prime[0], &inc);

    typename orbital_configuration_t::iterator it;
    it=segments_old.begin();
    for (int i=0; i<M.get_current_size(); i++) {
      R[i] = HYBRIDIZATION_TOOLS::interpolate_F(coor, new_segment.t_end()-it->t_start(), F);
      it++;
    }

    double det_rat = BLAS::ddot_(&M.get_current_size(),&R[0],&inc,&Q_prime[0],&inc);

    {
      overlap = 1;
      // take care of sign changes produced by segments which "wind around"
      if (k==(int)segments_old.size()-1) {
	it--;
	// check if last segment has been shifted across beta
	if ((new_segment.t_end()-new_segment.t_start())*(it->t_end()-it->t_start())<0) {
	  det_rat *= -1;
	  overlap = -1;	  
	}
      }
	
      if (det_rat < 0) {
	det_rat_sign = -1;
	det_rat *= -1;
      }
      else 
	det_rat_sign = 1;
    }

    return det_rat;
  }

  /*! 
   * \ingroup  HYBRIDIZATION
   *
   * \brief    Calculates the new hybridization matrix for shifting a vertex end point using sherman-morrison equations (A.4-9).
   *
   */
  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename vertex_vertex_matrix_t, typename hybridization_function_t>
  void shift_segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::compute_M_shift(int this_flavor,
												  Hybridization_vertex& new_segment, 
												  int k, 
												  vertex_vertex_matrix_t& M,
												  orbital_configuration_t& segments_old, 
												  hybridization_function_t& F,
												  std::vector<double>&      R, 
												  std::vector<double>&      Q_prime, 
												  double det_rat)
  {
    static int inc = 1;
    static std::vector<double> R_prime (1, 0.);

    double S_prime = 1./det_rat;

    if(M.get_current_size()>(int)R_prime.size())
      R_prime.resize(M.get_current_size());

    memset(&M(0,k), 0, sizeof(double)*M.get_current_size());
  
    HYBRIDIZATION_TOOLS::compute_R_prime(R, M, R_prime);
    R[k] = 0;

    HYBRIDIZATION_TOOLS::compute_M(Q_prime, R_prime, S_prime, M);

    BLAS::daxpy_(&M.get_current_size(), &S_prime, &Q_prime[0], &inc, &M(0,k), &inc);
  }  

}

#endif
