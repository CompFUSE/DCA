//-*-C++-*-

#ifndef SS_HYBRIDIZATION_TOOLS_H
#define SS_HYBRIDIZATION_TOOLS_H

namespace DCA
{  
/*!
 *  \class   HYBRIDIZATION_TOOLS
 *  \ingroup HYBRIDIZATION
 *
 *  \author  Bart Ydens
 *  \brief   This class implements the help functions for the insertion and removal of (anti-)segments. The help functions
 *  include the calculation of the deteminant ratio and the computation of the new hybridization matrix, using sherman-morrison equations.
 *
 */
  struct HYBRIDIZATION_TOOLS
  {
#include "type_definitions.h"

  public :

    static int cycle(int i, int size);

    template<typename Hybridization_function_t> 
    static double interpolate_F(int* coor_flavor, double tau, Hybridization_function_t& F);

    static double compute_length(double r, double l_max, double mu);

    template<class S> 
    static double compute_overlap(Hybridization_vertex segment, S& other_segments, int other_full_line, double BETA);
  
    template<class S> 
    static double segment_overlap(Hybridization_vertex segment, S& other_segments, int other_full_line, double BETA); 

    template <typename orbital_configuration_t> 
    static void compute_intervals(double t, 
				  double BETA, 
				  double& t_up, 
				  double& t_down, 
				  orbital_configuration_t& segments, 
				  typename orbital_configuration_t::iterator& s_up, 
				  typename orbital_configuration_t::iterator& s_down); 

    template <typename G, typename vertex_vertex_matrix_type, typename orbital_configuration_t> 
    static double det_rat_up(int                        this_flavor,
			     Hybridization_vertex&      new_segment, 
			     vertex_vertex_matrix_type& M, 
			     orbital_configuration_t&   segments_old, 
			     G& F, 
			     std::vector<double>& R, 
			     std::vector<double>& Q,
			     double& det_rat_sign, 
			     double& overlap);

    template <class G, typename vertex_vertex_matrix_type, typename orbital_configuration_t> 
    static void compute_M_up(int r, int s,
			     vertex_vertex_matrix_type& M, 
			     orbital_configuration_t&   segments_old, 
			     G& F, 
			     std::vector<double>& Fs, 
			     std::vector<double>& Fe,
			     double det_rat); 

    template <typename vertex_vertex_matrix_type, typename orbital_configuration_t> 
    static double det_rat_down(int r, int s, 
			       vertex_vertex_matrix_type& M, 
			       orbital_configuration_t& segments_old, 
			       double & det_rat_sign); 

    template<typename vertex_vertex_matrix_type> 
    static void compute_M_down(int r, int s, vertex_vertex_matrix_type& M);

    template <typename G, typename vertex_vertex_matrix_type, typename orbital_configuration_t> 
    static double det_rat_shift(int                        this_flavor,
				Hybridization_vertex&      new_segment,    
				int k,
				vertex_vertex_matrix_type& M, 
				orbital_configuration_t&   segments_old, 
				G& F, 
				std::vector<double>& R, 
				std::vector<double>& Q,
				double& det_rat_sign, 
				double& overlap);

    template <typename vertex_vertex_matrix_type> 
    static void compute_M_shift(int k, 
				vertex_vertex_matrix_type& M,
				std::vector<double>& Fs, 
				std::vector<double>& Fe,
				double det_rat); 

    template<typename vertex_vertex_matrix_type>
    static void compute_Q_prime(std::vector<double>& Q, vertex_vertex_matrix_type& M, std::vector<double>& Q_prime);
    template<typename vertex_vertex_matrix_type>
    static void compute_R_prime(std::vector<double>& R, vertex_vertex_matrix_type& M, std::vector<double>& R_prime);
    template <typename vertex_vertex_matrix_type> 
    static void compute_M(std::vector<double>& Q_prime, std::vector<double>& R_prime, double S_prime, vertex_vertex_matrix_type& M);

  };

  int HYBRIDIZATION_TOOLS::cycle(int i, int size) 
  {
    return (i>0 ? i-1 : size-1); 
  }

  template<typename Hybridization_function_t> 
  double HYBRIDIZATION_TOOLS::interpolate_F(int* coor_flavor, double tau, Hybridization_function_t& F) 
  {
    static int* coor = new int[F.signature()];

    coor[0] = coor_flavor[0];
    coor[1] = coor_flavor[1];
    coor[2] = coor_flavor[0];
    coor[3] = coor_flavor[1];
    coor[4] = 0;

    double sign=1;

    if (tau<0) {
      tau += time_domain_type::get_beta();
      sign=-1;
    }
    
    return sign*polynomial_interpolation<t,1>::execute(tau, F, coor);
  }

  double HYBRIDIZATION_TOOLS::compute_length(double r, double l_max, double mu) 
  {
    if (mu == 0)
      return r*l_max;	
    else
      return 1/mu*log(r*(exp(mu*l_max)-1)+1);
  }

  template <class S> 
  double HYBRIDIZATION_TOOLS::compute_overlap(Hybridization_vertex segment, S& other_segments, int other_full_line, double BETA) 
  {
    if (segment.t_start()<segment.t_end())
      return segment_overlap(segment, other_segments, other_full_line, BETA);
    else
      {
	double other_length=0;
	Hybridization_vertex segment1(0,segment.t_end());
	Hybridization_vertex segment2(segment.t_start(), BETA);
	other_length += segment_overlap(segment1, other_segments, other_full_line, BETA);
	other_length += segment_overlap(segment2, other_segments, other_full_line, BETA);
	return other_length;
      }
  }

  template <class S> 
  double HYBRIDIZATION_TOOLS::segment_overlap(Hybridization_vertex segment, S& other_segments, int other_full_line, double BETA) 
  {
    double length = (segment.t_start()<segment.t_end() ? segment.t_end()-segment.t_start() : segment.t_end()-segment.t_start()+BETA);
    double t_final = segment.t_start()+length;
    double t = segment.t_start();
    double t_final_segment;		
    double other_length=0;
    if (other_full_line==1)
      other_length=length;
    else if (other_segments.size()>0){
      typename S::iterator it;
      it = lower_bound(other_segments.begin(), other_segments.end(), t);	
		
      if (it!=other_segments.begin()) {
	it--;
	t_final_segment = (it->t_start()<it->t_end() ? it->t_end() : it->t_end()+BETA);
	if (t<t_final_segment) {
	  other_length += (t_final_segment<t_final ? t_final_segment-t : t_final-t);
	}
	it++;

      }
      while(it!=other_segments.end() && it->t_start()<t_final) {
	t_final_segment = (it->t_start()<it->t_end() ? it->t_end() : it->t_end()+BETA);
	other_length += (t_final_segment<t_final ? t_final_segment-it->t_start() : t_final-it->t_start());
	it++;
      }
      // check if last segment overlaps
      it=other_segments.end();
      it--;
      if (it->t_end()<it->t_start() && t<it->t_end()) {
	other_length += (t_final<it->t_end() ? t_final-t : it->t_end()-t);
      }
    } 	
    return other_length;   
  }

  template <typename orbital_configuration_t> 
  void HYBRIDIZATION_TOOLS::compute_intervals(double t, 
					      double BETA, 
					      double& t_up, 
					      double& t_down, 
					      orbital_configuration_t& segments, 
					      typename orbital_configuration_t::iterator& s_up, 
					      typename orbital_configuration_t::iterator& s_down) 
  {  
    if (segments.size() == 0) 
      {
	t_up = BETA;
	t_down = BETA;
	s_up = segments.end();
	s_down = segments.end();
      }
    else 
      {
	s_up = lower_bound(segments.begin(), segments.end(), t);
	
	if (s_up == segments.begin()) {
	  s_down = segments.end(); s_down--;
	  if (s_down->t_end() < s_down->t_start())
	    t_down = t - s_down->t_end();
	  else
	    t_down = t + BETA - s_down->t_end();
	}
	else {
	  s_down = s_up; s_down--;
	  if (s_down->t_end()>s_down->t_start())
	    t_down = t - s_down->t_end();
	  else 
	    t_down = t - (BETA+s_down->t_end());
	}
	
	if(s_up == segments.end()) {
	  t_up = BETA - t + segments.begin()->t_start();
	}
	else {
	  t_up = s_up->t_start() - t;
	}
    }
  }


  /*! 
   * \ingroups  HYBRIDIZATION TOOLS
   *
   * \brief    Calculates the determinant ratio for inserting a new vertex. The determinant ratio is given by (A.10)
   * \f{eqnarray*}{
   * \frac{det(M^{k+1}_{\sigma})}{det(M^{k}_{\sigma})} = det(S - R M^k Q)
   * \f}
   * with\f$ S = F(\tau^n_{e} - \tau^n_{s})\f$, 
   * R a (1 x k)-vector with \f$R[i] =  F(\tau^n_{e} - \tau^i_{s})\f$,
   * Q a (k x 1)-vector with \f$Q[i] =  F(\tau^i_{e} - \tau^n_{s})\f$. 
   *
   */
  template <typename Hybridization_function_t, typename vertex_vertex_matrix_type, typename orbital_configuration_t> 
  double HYBRIDIZATION_TOOLS::det_rat_up(int                        this_flavor,
					 Hybridization_vertex&      new_segment, 
					 vertex_vertex_matrix_type& M, 
					 orbital_configuration_t&   segments_old, 
					 Hybridization_function_t&  F, 
					 std::vector<double>&       R, 
					 std::vector<double>&       Q_prime,
					 double&                    det_rat_sign, 
					 double&                    overlap) 
  {
    static int inc = 1;
    static std::vector<double> Q      (1, 0.);
    if(M.get_current_size()>(int)Q.size())
      Q.resize(M.get_current_size());

    static int* coor = new int[2];
    static nu nu_obj;
    nu_obj.linind_2_subind(this_flavor, coor);

    // S = F(\tau^n_{e} - \tau^n_{s})
    // R[i] =  F(\tau^n_{e} - \tau^i_{s})
    // Q[i] =  F(\tau^i_{e} - \tau^n_{s})
    double det_rat = interpolate_F(coor, new_segment.t_end()-new_segment.t_start(), F);

    if(M.get_current_size() > 0){
      typename orbital_configuration_t::iterator it=segments_old.begin();
      for (size_t i=0; i<segments_old.size(); i++) 
	{
	  R[i] = interpolate_F(coor, new_segment.t_end()-it->t_start(), F);
	  Q[i] = interpolate_F(coor, it->t_end()-new_segment.t_start(), F);
	  it++;
	}

      compute_Q_prime(Q, M, Q_prime);

      int size = M.get_current_size();
      det_rat += BLAS::ddot_(&size, &R[0], &inc, &Q_prime[0], &inc);
    }
    
    // take care of sign changes produced by segments which "wind around"
    {
      if(new_segment.t_end() < new_segment.t_start()) {
	det_rat *= -1;	  
	overlap  = -1;	
      }
      else
	overlap  = 1;
      
      if (det_rat < 0) {
	det_rat     *= -1;
	det_rat_sign = -1;
      }
      else 
	det_rat_sign = 1;
    }

    return det_rat;
  }


  /*! 
   * \ingroups  HYBRIDIZATION TOOLS
   *
   * \brief    Calculates the new hybridization matrix for inserting a new vertex using sherman-morrison equations (A.4-9).
   *
   */
  template <typename Hybridization_function_t, 
	    typename vertex_vertex_matrix_type, 
	    typename orbital_configuration_t> 
  void HYBRIDIZATION_TOOLS::compute_M_up(int r, int s,
					 vertex_vertex_matrix_type& M, 
					 orbital_configuration_t&   segments_old, 
					 Hybridization_function_t&  F, 
					 std::vector<double>& R, 
					 std::vector<double>& Q_prime,
					 double S_prime_inv) 
  {
    static std::vector<double> R_prime (1, 0.);

    int i_new, j_new;
    int size=M.get_current_size();

    if(size > 0){
      if(size > (int)R_prime.size())
	R_prime.resize(M.get_current_size());

      compute_R_prime(R, M, R_prime);

      compute_M(Q_prime, R_prime, 1./S_prime_inv, M);
    }

    if(r==0 && s!=0) // need to permute indices of R, L, M
      M.cycle_column_forward();	

    if(r < M.get_current_size() && s < M.get_current_size())
      M.insert_row_and_column(r,s);
    else if(r < M.get_current_size() && s == M.get_current_size())
      M.insert_row_and_add_column(r);
    else if(r == M.get_current_size() && s < M.get_current_size())
      M.add_row_and_insert_column(s);
    else{
      assert(r == s && r == size);
      M.resize(size+1);
    }

    // row k+1 and column k
    if (r!=0 || r == s) { // segments remain in the usual order
      for (int i=0; i<size; i++) {
	i_new = (i<r ? i : i+1);
	j_new = (i<s ? i : i+1);
	
	M(i_new,s) = Q_prime[i]/S_prime_inv;
	M(r,j_new) = R_prime[i]/S_prime_inv;
      }
    }
    else { // need to permute indices of R, L, M
      for (int i=0; i<size; i++) {
	i_new = (i<r ? i : i+1);
	j_new = (i<s ? i : i+1);
	
	M(i_new,s) = Q_prime[i]/S_prime_inv;
	M(r,j_new) = R_prime[HYBRIDIZATION_TOOLS::cycle(i,size)]/S_prime_inv;
      }
    }

    // fill S_prime
    M(r,s) = 1./S_prime_inv;
  }  

  /*! 
   * \ingroups  HYBRIDIZATION TOOLS
   *
   * \brief    Calculates the determinant ratio for removing a vertex. 
   *
   */
  template <typename vertex_vertex_matrix_type, typename orbital_configuration_t> 
  double HYBRIDIZATION_TOOLS::det_rat_down(int r, int s, 
					   vertex_vertex_matrix_type& M, 
					   orbital_configuration_t& segments_old, 
					   double & det_rat_sign) 
  {
    double det_rat = M(r,s);

    // take care of sign changes produced by segments which "wind around"
    if (r==int(segments_old.size())-1) {
      typename orbital_configuration_t::iterator it=segments_old.end(); it--;
      if (it->t_end() < it->t_start())
	det_rat *= -1;	  
    }
	
    if (det_rat < 0) {
      det_rat_sign = -1;
      det_rat *= -1;
    }
    else {
      det_rat_sign = 1;
    }
  
    return det_rat;
  }

  /*! 
   * \ingroups  HYBRIDIZATION TOOLS
   *
   * \brief    Calculates the new hybridization matrix for removing a vertex using sherman-morrison equations (A.4-9).
   *
   */
  template<typename vertex_vertex_matrix_type>
  void HYBRIDIZATION_TOOLS::compute_M_down(int r, int s, vertex_vertex_matrix_type& M) 
  {
    static int inc = 1;
    
    static std::vector<double> Q_prime (1, 0.);
    static std::vector<double> R_prime (1, 0.);

    if(M.get_current_size() > 1)
      {
	if(M.get_current_size()>(int)R_prime.size()){
	  Q_prime.resize(M.get_current_size());
	  R_prime.resize(M.get_current_size());
	}

	int incy = M.get_global_size();
	int size = M.get_current_size();
	BLAS::dcopy_(&size, &M(0,s), &inc , &Q_prime[0], &inc);
	BLAS::dcopy_(&size, &M(r,0), &incy, &R_prime[0], &inc);

	compute_M(Q_prime,R_prime,-1./Q_prime[r],M);
	M.remove_row_and_column(r,s);
	
	if(r==0 && s != 0) // need to permute indices of R, L, M
	  M.cycle_column_backward();
      }
    else
      M.resize(0);
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
  template <typename Hybridization_function_t, typename vertex_vertex_matrix_type, typename orbital_configuration_t> 
  double HYBRIDIZATION_TOOLS::det_rat_shift(int this_flavor,
					    Hybridization_vertex& new_segment, 
					    int k, 
					    vertex_vertex_matrix_type& M,
					    orbital_configuration_t& segments_old, 
					    Hybridization_function_t& F, 
					    std::vector<double>&      R, 
					    std::vector<double>&      Q_prime,
					    double & det_rat_sign, 
					    double & overlap) 
  {
    static int inc   = 1;
    static int* coor = new int[2];
    static nu nu_obj;
    nu_obj.linind_2_subind(this_flavor, coor);

    int size = M.get_current_size();

    BLAS::dcopy_(&size, &M(0,k), &inc, &Q_prime[0], &inc);

    typename orbital_configuration_t::iterator it;
    it=segments_old.begin();
    for (int i=0; i<M.get_current_size(); i++) {
      R[i] = interpolate_F(coor, new_segment.t_end()-it->t_start(), F);
      it++;
    }

    double det_rat = BLAS::ddot_(&size,&R[0],&inc,&Q_prime[0],&inc);

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
  template <typename vertex_vertex_matrix_type>
  void HYBRIDIZATION_TOOLS::compute_M_shift(int k, 
					    vertex_vertex_matrix_type& M,
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
  
    compute_R_prime(R, M, R_prime);
    R[k] = 0;

    compute_M(Q_prime, R_prime, S_prime, M);

    int size = M.get_current_size();

    BLAS::daxpy_(&size, &S_prime, &Q_prime[0], &inc, &M(0,k), &inc);
  }  

  /*! 
   * \ingroups  HYBRIDIZATION TOOLS
   *
   * \brief    Calculates Q' = -M*Q
   *
   */
  template<typename vertex_vertex_matrix_type>
  void HYBRIDIZATION_TOOLS::compute_Q_prime(std::vector<double>& Q, vertex_vertex_matrix_type& M,std::vector<double>& Q_prime) 
  {
    static gemv_plan  <double> gemv_pln(1, 1);
    {
      gemv_pln.TRANS         = 'N';
      gemv_pln.alpha         = -1;
      gemv_pln.M             = M.get_current_size();
      gemv_pln.N             = M.get_current_size();
      gemv_pln.LDA           = M.get_global_size();
      gemv_pln.matrix        = &M(0,0);
      gemv_pln.vector_source = &Q[0];
      gemv_pln.vector_target = &Q_prime[0];

      gemv_pln.execute_plan();
    }
  }

  /*! 
   * \ingroups  HYBRIDIZATION TOOLS
   *
   * \brief     Calculates R' = -R*M
   *
   */
  template<typename vertex_vertex_matrix_type>
  void HYBRIDIZATION_TOOLS::compute_R_prime(std::vector<double>& R, vertex_vertex_matrix_type& M,std::vector<double>& R_prime) 
  {
    static gemv_plan  <double> gemv_pln(1, 1);
    
    gemv_pln.TRANS         = 'T';
    gemv_pln.alpha         = -1;
    gemv_pln.M             = M.get_current_size();
    gemv_pln.N             = M.get_current_size();
    gemv_pln.LDA           = M.get_global_size();
    gemv_pln.matrix        = &M(0,0);
    gemv_pln.vector_source = &R[0];
    gemv_pln.vector_target = &R_prime[0];
    
    gemv_pln.execute_plan();
  }


  /*! 
   * \ingroups  HYBRIDIZATION TOOLS
   *
   * \brief    Calculates new M matrix M_new = Q'*R'*S'
   *
   */
  template<typename vertex_vertex_matrix_type>
  void HYBRIDIZATION_TOOLS::compute_M(std::vector<double>& Q_prime,std::vector<double>& R_prime, double S_prime, vertex_vertex_matrix_type& M) 
  {
    static gemm_plan  <double> gemm_pln;

    gemm_pln.alpha         = S_prime;
    gemm_pln.beta          = 1;

    gemm_pln.M             = M.get_current_size();
    gemm_pln.N             = M.get_current_size();
    gemm_pln.K             = 1;

    gemm_pln.LDA           = M.get_global_size();
    gemm_pln.LDB           = 1;
    gemm_pln.LDC           = M.get_global_size();

    gemm_pln.A             = &Q_prime[0];
    gemm_pln.B             = &R_prime[0];
    gemm_pln.C             = &M(0,0);

    gemm_pln.execute_plan();

  }


}

#endif
