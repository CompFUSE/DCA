//-*-C++-*-

#ifndef SEGMENT_TOOLS_H
#define SEGMENT_TOOLS_H

namespace QMC {
  
/*!
 *  \ingroup HYBRIDIZATION
 *
 *  \author  Bart Ydens
 *  \brief   This class implements the insertion and removal of segments \f$(c^{\dagger}-c\ pair)\f$.
 *
 *  \f{eqnarray}{
 *   W_{acc}(c_k \rightarrow c_{k+1}) = min \left(1, \frac{\beta l_{max}}{k+1} \frac{det F(k+1)}{det F(k)}\frac{W_{loc}(k+1)}{W_{loc}(k)} \right)
 *  \f}
 *
 */
  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  class segment_tools
  {
#include "type_definitions.h"
   
    typedef typename configuration_t::orbital_configuration_t orbital_configuration_t;

    typedef segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t> THIS_TYPE;

  public:

    segment_tools(configuration_t& configuration,
		  parameters_t&    parameters,
		  MOMS_t&          MOMS,
		  concurrency_t&   concurrency);

    ~segment_tools();

    template<typename function_type_1, typename function_type_2>
    double insert_segment(int this_flavor, 
			double mu,
			double& sign,
			function_type_1& M,
			function_type_2& F);
    
    template<typename function_type_1, typename function_type_2>
    double remove_segment(int this_flavor, 
			double mu,
			double& sign,
			function_type_1& M,
			function_type_2& F);

  private:
    configuration_t& configuration;
    
    concurrency_t&   concurrency;
    parameters_t&    parameters;
    MOMS_t&          MOMS;

    int    FLAVORS;
    double BETA;
  };

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::segment_tools(configuration_t& configuration_ref,
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
  segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::~segment_tools()
  {}

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename function_type_1, typename function_type_2>
  double segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::insert_segment(int this_flavor, 
											   double mu,
											   double& sign,
											   function_type_1& M,
											   function_type_2& F)
  {
//     cout << __FUNCTION__ << endl;
    assert(this_flavor >=0 && this_flavor < FLAVORS);

    double t = BETA*concurrency.get_random_number();

    double t_up; // distance to next segment up
    double t_down; // distance to next segment down

    orbital_configuration_t& vertices = configuration.get_vertices(this_flavor);

    typename orbital_configuration_t::iterator s_up; // iterator of the segment up
    typename orbital_configuration_t::iterator s_down; // iterator of the segment down

    HYBRIDIZATION_TOOLS::compute_intervals(t, BETA, t_up, t_down, vertices, s_up, s_down);
    
    double trace = -1;

    if (t_down>0) 
      { // t does not lie on a segment -> it's possible to insert a new one starting from t
	double length = HYBRIDIZATION_TOOLS::compute_length(concurrency.get_random_number(), t_up, 0);
	  
	Hybridization_vertex segment_insert;
	segment_insert.set_t_start(t);
	double t_final = t + length;
	if (t_final > BETA)
	  segment_insert.set_t_end(t_final-BETA);
	else
	  segment_insert.set_t_end(t_final);
	
	double otherlength_u=0;

	for (int i=0; i<FLAVORS; i++) {
	  if (i==this_flavor) continue;
	  double other_length = HYBRIDIZATION_TOOLS::compute_overlap(segment_insert, configuration.get_vertices(i), configuration.get_full_line(i), BETA);//compute_overlap(segment_insert, other_segments[i], other_full_line[i], BETA);
	  otherlength_u += other_length*MOMS.H_interactions(i, this_flavor, 0);//u(i, this_flavor);
	}
	  
	double log_prob, overlap, det_rat, det_rat_sign;
	std::vector<double> R(vertices.size()), Q(vertices.size());
	  
	det_rat = HYBRIDIZATION_TOOLS::det_rat_up(this_flavor, segment_insert, M(this_flavor), vertices, F, R, Q, det_rat_sign, overlap);

	log_prob = log(BETA*t_up/(vertices.size()+1)*det_rat)+mu*length-otherlength_u;
	  
	if (log(concurrency.get_random_number()) < log_prob) {

	  int position=0;
	  for (typename orbital_configuration_t::iterator it=vertices.begin(); it!=s_up; it++)
	    position++;

	  HYBRIDIZATION_TOOLS::compute_M_up(segment_insert, position, M(this_flavor), vertices, F, R, Q, det_rat*overlap);
	  sign *= det_rat_sign;
	  vertices.insert(s_up, segment_insert);
	  trace = std::exp(mu*length-otherlength_u);
	}
      }
    return trace;
  }

  template<typename configuration_t, typename parameters_t, typename MOMS_t, typename concurrency_t>
  template<typename function_type_1, typename function_type_2>
  double segment_tools<configuration_t, parameters_t, MOMS_t, concurrency_t>::remove_segment(int this_flavor, 
											   double mu,
											   double& sign,
											   function_type_1& M,
											   function_type_2& F)
  {
//     cout << __FUNCTION__ << endl;
    assert(this_flavor >=0 && this_flavor < FLAVORS);

    orbital_configuration_t& vertices = configuration.get_vertices(this_flavor);

    if (vertices.size()==0)
      return -1;

    typename orbital_configuration_t::iterator s_up;   // iterator of the segment up
    typename orbital_configuration_t::iterator s_down; // iterator of the segment down

    int position = concurrency.get_random_number()*vertices.size();
    s_down = vertices.begin();
    for (int i=0; i<position; i++)
      s_down++;
    s_up=s_down;
    s_up++;
    if (s_up==vertices.end())
      s_up = vertices.begin();

    double length = s_down->t_end()-s_down->t_start();
    if (length < 0) length += BETA;
	  
    double t_total = s_up->t_start()-s_down->t_start();
    if (t_total <= 0) t_total += BETA;
      
    Hybridization_vertex segment_remove = *s_down;
	
    double otherlength_u=0;
    for (int i=0; i<FLAVORS; i++) {
      if (i==this_flavor) continue;
      double other_length = HYBRIDIZATION_TOOLS::compute_overlap(segment_remove, configuration.get_vertices(i), configuration.get_full_line(i), BETA);//compute_overlap(segment_remove, other_segments[i], other_full_line[i], BETA);
      otherlength_u += other_length*MOMS.H_interactions(i, this_flavor, 0);//u(i, this_flavor);
    }

    double log_prob/*, overlap*/, det_rat, det_rat_sign;
	
    det_rat = HYBRIDIZATION_TOOLS::det_rat_down(position, M(this_flavor), vertices, det_rat_sign);	  
	
    log_prob = log(BETA*t_total/vertices.size()/det_rat)+length*mu-otherlength_u;
	  
    if (log(concurrency.get_random_number()) < -log_prob) {
      HYBRIDIZATION_TOOLS::compute_M_down(position, M(this_flavor));	
      sign *= det_rat_sign;
      vertices.erase(s_down);
      return std::exp(-mu*length+otherlength_u);
    }
    return -1;
  }
}

#endif
