//-*-C++-*-

#ifndef DCA_CLUSTER_EXCLUSION_STEP_H
#define DCA_CLUSTER_EXCLUSION_STEP_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace DCA
{
  template<typename parameters_type, typename MOMS_type>
  class cluster_exclusion
  {

    typedef typename parameters_type::profiler_type    profiler_t;
    typedef typename parameters_type::concurrency_type concurrency_type;

  public:

    cluster_exclusion(parameters_type&    parameters_ref,
		      MOMS_type& MOMS_ref);

    ~cluster_exclusion();
  
    void execute();

  private:

    void compute_G0_K_w_cluster_excluded();
    void compute_G0_R_t_cluster_excluded();

    void plot_G0_R_t_cluster_excluded();

  private:

    parameters_type& parameters;
    MOMS_type&       MOMS;
  };

  template<typename parameters_type, typename MOMS_type>
  cluster_exclusion<parameters_type, MOMS_type>::cluster_exclusion(parameters_type& parameters_ref,
								   MOMS_type&       MOMS_ref):
    parameters(parameters_ref),
    MOMS(MOMS_ref)
  {}

  template<typename parameters_type, typename MOMS_type>
  cluster_exclusion<parameters_type, MOMS_type>::~cluster_exclusion()
  {}

  template<typename parameters_type, typename MOMS_type>
  void cluster_exclusion<parameters_type, MOMS_type>::execute()
  {
    compute_G0_K_w_cluster_excluded();
  
    compute_G0_R_t_cluster_excluded();
  }

  /*
   *   G = G_0 + G_0*S*G = G_0 * (1 + S*G)
   *
   *   G_0 = G*(1 + S*G)^-1
   */
  template<typename parameters_type, typename MOMS_type>
  void cluster_exclusion<parameters_type, MOMS_type>::compute_G0_K_w_cluster_excluded()
  {
    profiler_t profiler(__FUNCTION__, "cluster_exclusion", __LINE__);

    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> G_matrix("G_matrix" , nu::dmn_size());
    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> S_matrix("S_matrix", nu::dmn_size());
    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> one_plus_S_G("one_plus_S_G", nu::dmn_size());
    LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> G0_matrix("G0_matrix", nu::dmn_size());

    LIN_ALG::GEINV<LIN_ALG::CPU>::plan<std::complex<double> > geinv_obj(one_plus_S_G);

    for(int w_ind=0; w_ind<w::dmn_size(); w_ind++){
      for(int K_ind=0; K_ind<k_DCA::dmn_size(); K_ind++){
	
	for(int j=0; j<nu::dmn_size(); j++)
	  for(int i=0; i<nu::dmn_size(); i++)
	    G_matrix(i,j) = MOMS.G_k_w(i,j,K_ind,w_ind);
      
	for(int j=0; j<nu::dmn_size(); j++)
	  for(int i=0; i<nu::dmn_size(); i++)
	    S_matrix(i,j) = MOMS.Sigma_cluster(i,j,K_ind,w_ind);

	for(int j=0; j<nu::dmn_size(); j++)
	  for(int i=0; i<nu::dmn_size(); i++)
	    one_plus_S_G(i,j) = 0;

	for(int j=0; j<nu::dmn_size(); j++)
	  for(int i=0; i<nu::dmn_size(); i++)
	    for(int l=0; l<nu::dmn_size(); l++)
	      one_plus_S_G(i,j) += S_matrix(i,l)*G_matrix(l,j);

	for(int i=0; i<nu::dmn_size(); i++)
	  one_plus_S_G(i,i) += 1.;

	geinv_obj.execute(one_plus_S_G);

	for(int j=0; j<nu::dmn_size(); j++)
	  for(int i=0; i<nu::dmn_size(); i++)
	    G0_matrix(i,j) = 0;

	for(int j=0; j<nu::dmn_size(); j++)
	  for(int i=0; i<nu::dmn_size(); i++)
	    for(int l=0; l<nu::dmn_size(); l++)
	      G0_matrix(i,j) += G_matrix(i,l)*one_plus_S_G(l,j);

	for(int j=0; j<nu::dmn_size(); j++)
	  for(int i=0; i<nu::dmn_size(); i++)
	    MOMS.G0_k_w_cluster_excluded(i,j,K_ind,w_ind) = G0_matrix(i,j);
      }
    }
  }

  template<typename parameters_type, typename MOMS_type>
  void cluster_exclusion<parameters_type, MOMS_type>::compute_G0_R_t_cluster_excluded()
  {
    profiler_t profiler(__FUNCTION__, "cluster_exclusion", __LINE__);

    MOMS.G0_k_w_cluster_excluded -= MOMS.G0_k_w;

    {
      math_algorithms::functional_transforms::TRANSFORM<w, t>::execute(MOMS.G0_k_w_cluster_excluded, 
						MOMS.G0_k_t_cluster_excluded);

      MOMS.G0_k_t_cluster_excluded += MOMS.G0_k_t;

      math_algorithms::functional_transforms::TRANSFORM<k_DCA, r_DCA>::execute(MOMS.G0_k_t_cluster_excluded, 
							MOMS.G0_r_t_cluster_excluded);
    }

    MOMS.G0_k_w_cluster_excluded += MOMS.G0_k_w;
  }

  template<typename parameters_type, typename MOMS_type>
  void cluster_exclusion<parameters_type, MOMS_type>::plot_G0_R_t_cluster_excluded()
  {
    {
      FUNC_LIB::function<float, t> tmp("G0_k_t");
    
      Gnuplot plot_obj("lines");
      for(int R_ind=0; R_ind<r_DCA::dmn_size(); R_ind++)  
	{
	  for(int t_ind=0; t_ind<t::dmn_size(); t_ind++)
	    tmp(t_ind) = MOMS.G0_k_t(0,0,R_ind,t_ind);
	
	  SHOW::execute(plot_obj, tmp);      
	}
    
      plot_obj.showonscreen();
    }

    {
      FUNC_LIB::function<float, t> tmp("G0_k_t_cluster_excluded");
    
      Gnuplot plot_obj("lines");
      for(int R_ind=0; R_ind<r_DCA::dmn_size(); R_ind++)  
	{
	  for(int t_ind=0; t_ind<t::dmn_size(); t_ind++)
	    tmp(t_ind) = MOMS.G0_k_t_cluster_excluded(0,0,R_ind,t_ind);
	
	  SHOW::execute(plot_obj, tmp);      
	}
    
      plot_obj.showonscreen();
    }

    {
      FUNC_LIB::function<float, t> tmp("G0_r_t");
    
      Gnuplot plot_obj("lines");
      for(int R_ind=0; R_ind<r_DCA::dmn_size(); R_ind++)  
	{
	  for(int t_ind=0; t_ind<t::dmn_size(); t_ind++)
	    tmp(t_ind) = MOMS.G0_r_t_cluster_excluded(0,0,R_ind,t_ind);
	
	  SHOW::execute(plot_obj, tmp);      
	}
    
      plot_obj.showonscreen();
    }
  }

}

#endif
