//-*-C++-*-

#ifndef DCA_COARSEGRAINING_ROUTINES_H
#define DCA_COARSEGRAINING_ROUTINES_H

namespace DCA
{

  template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  struct coarsegraining_functions
  {
#include "type_definitions.h"

    int K_ind;
    int w_ind;

    scalar_type mu;

    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >* H_k;
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >* A_k;

    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> > I_q;
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> > H_q;
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> > A_q;
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> > S_q;
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> > G_q;
  };


  template<typename parameters_type, typename K_dmn>
  class coarsegraining_routines
  {
#include "type_definitions.h"

    typedef typename K_dmn::parameter_type k_cluster_type;

    const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

    typedef typename parameters_type::concurrency_type concurrency_type;

    typedef dmn_0<MATH_ALGORITHMS::tetrahedron_mesh<K_dmn> >             tetrahedron_dmn;

    typedef MATH_ALGORITHMS::gaussian_quadrature_domain<tetrahedron_dmn> quadrature_dmn;

  public:

    coarsegraining_routines(parameters_type& parameters_ref);

    ~coarsegraining_routines();

  protected:

    void compute_tetrahedron_mesh(int k_mesh_refinement,
                                  int number_of_periods);

    void compute_gaussian_mesh(int k_mesh_refinement,
                               int gaussian_quadrature_rule,
                               int number_of_periods);

    template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
    void wannier_interpolation(int K_ind,
                               FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& f_k,
                               FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& f_q);

    template<typename scalar_type, typename tmp_scalar_type, typename q_dmn_t>
    void compute_I_q(tmp_scalar_type value,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q);

    template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
    void compute_H_q(int K_ind,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& H_k,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q);

    /***********************************************************************************
     ***
     ***         Routines for DCA
     ***
     ***********************************************************************************/

    template<typename scalar_type, typename q_dmn_t>
    void compute_S_q(int K_ind, int w_ind,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn  , w> >& S_K_w,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t   > >& S_q);

    template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
    void compute_G_q_w(int K_ind, int w_ind,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t   > >& H_k,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn  , w> >& S_K,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t   > >& I_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t   > >& H_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t   > >& S_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t   > >& G_q);

    template<typename scalar_type, typename q_dmn_t>
    void compute_S_q(int K_ind, int w_ind,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn  , w_REAL> >& S_K_w,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t   > >& S_q);

    template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
    void compute_G_q_w(int K_ind, int w_ind,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t        > >& H_k,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn  , w_REAL> >& S_K,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t        > >& I_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t        > >& H_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t        > >& S_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t        > >& G_q);

    template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
    void compute_G_q_t(int K_ind, int t_ind,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& H_k,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q);

    /***********************************************************************************
     ***                                                                             ***
     ***         Routines for DCA+
     ***                                                                             ***
     ***********************************************************************************/

    template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
    void compute_S_q_from_A_k(int K_ind, int w_ind,
                              FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& A_k,
                              FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& A_q,
                              FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& S_q);

    template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
    void compute_G_q_w(int K_ind, int w_ind,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& H_k,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& A_k,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& A_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& S_q,
                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q);

    template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
    void compute_G_q_w(coarsegraining_functions<scalar_type, k_dmn_t, q_dmn_t>& coarsegraining_functions_ref);

    /***********************************************************************************
     ***                                                                             ***
     ***         tetrahedron-integration
     ***                                                                             ***
     ***********************************************************************************/

    /*
    typedef dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_K     > > tet_dmn_type;
    typedef dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN> > tet_0_dmn_type;

    template<typename scalar_type>
    struct tetrahedron_integration_functions
    {
      int size;

      FUNC_LIB::function<             scalar_type , tet_dmn_type>*                 w_tet_ptr;
      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >* G_tet_ptr;

      std::vector<FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> > >  G_int_vec;
    };

    template<typename scalar_type, typename tet_dmn_t>
    void tetrahedron_integration_st(FUNC_LIB::function<scalar_type              , tet_dmn_t>&                 w_tet,
                                    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                    FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int);

    template<typename scalar_type, typename tet_dmn_t>
    void tetrahedron_integration_mt(FUNC_LIB::function<scalar_type              , tet_dmn_t>&                 w_tet,
                                    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                    FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int);

    template<typename scalar_type, typename tet_dmn_t>
    void tetrahedron_integration_st_1D(FUNC_LIB::function<scalar_type              , tet_dmn_t>&                 w_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int);

    template<typename scalar_type, typename tet_dmn_t>
    void tetrahedron_integration_mt_1D(FUNC_LIB::function<scalar_type              , tet_dmn_t>&                 w_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int);

    template<typename scalar_type, typename tet_dmn_t>
    void tetrahedron_integration_st_2D(FUNC_LIB::function<             scalar_type , tet_dmn_t>&                 w_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int);

    template<typename scalar_type>
    static void* threaded_tetrahedron_integration_2D(void* data);

    template<typename scalar_type, typename tet_dmn_t>
    void tetrahedron_integration_mt_2D(FUNC_LIB::function<scalar_type              , tet_dmn_t>&                 w_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int);

    template<typename scalar_type, typename tet_dmn_t>
    void tetrahedron_integration_st_3D(FUNC_LIB::function<scalar_type              , tet_dmn_t>&                 w_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int);

    template<typename scalar_type, typename tet_dmn_t>
    void tetrahedron_integration_mt_3D(FUNC_LIB::function<scalar_type              , tet_dmn_t>&                 w_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                       FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int);
    */

  protected:

    parameters_type&  parameters;
    concurrency_type& concurrency;
  };

  template<typename parameters_type, typename K_dmn>
  coarsegraining_routines<parameters_type, K_dmn>::coarsegraining_routines(parameters_type& parameters_ref):
    parameters (parameters_ref),
    concurrency(parameters.get_concurrency())
  {
  }

  template<typename parameters_type, typename K_dmn>
  coarsegraining_routines<parameters_type, K_dmn>::~coarsegraining_routines()
  {}

  template<typename parameters_type, typename K_dmn>
  void coarsegraining_routines<parameters_type, K_dmn>::compute_tetrahedron_mesh(int k_mesh_refinement,
                                                                                 int number_of_periods)
  {
    MATH_ALGORITHMS::tetrahedron_mesh<typename K_dmn::parameter_type> mesh(k_mesh_refinement);

    quadrature_dmn::translate_according_to_period(number_of_periods, mesh);

    {
      coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_size()     = 0;
      coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_weights() .resize(0);
      coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_elements().resize(0);

      for(int l=0; l<mesh.size(); l++)
        mesh[l].update_tetrahedron_domain(coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_size(),
                                          coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_weights(),
                                          coarsegraining_domain<K_dmn, TETRAHEDRON_K>::get_elements());


      coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_size()     = 0;
      coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_weights() .resize(0);
      coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_elements().resize(0);

      for(int l=0; l<mesh.size(); l++)
        mesh[l].update_tetrahedron_domain(coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_size(),
                                          coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_weights(),
                                          coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_elements());
    }

//     SHOW::plot_points(coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN>::get_elements());
//     throw std::logic_error(__FUNCTION__);
  }

  template<typename parameters_type, typename K_dmn>
  void coarsegraining_routines<parameters_type, K_dmn>::compute_gaussian_mesh(int k_mesh_refinement,
                                                                              int gaussian_quadrature_rule,
                                                                              int number_of_periods)
  {
    //     quadrature_dmn::initialize_Brillouin_zone(parameters.get_k_mesh_refinement(),
    //                                               parameters.get_gaussian_quadrature_rule(),
    //                                               parameters.get_number_of_periods());

    {
      quadrature_dmn::initialize_Brillouin_zone(k_mesh_refinement,
                                                gaussian_quadrature_rule,
                                                number_of_periods);

      //       SHOW::plot_points(K_dmn::get_elements());
      //       SHOW::plot_points(quadrature_dmn::get_elements());
    }

    {
      coarsegraining_domain<K_dmn, ORIGIN>::get_size()     = quadrature_dmn::get_size();
      coarsegraining_domain<K_dmn, ORIGIN>::get_weights()  = quadrature_dmn::get_weights();
      coarsegraining_domain<K_dmn, ORIGIN>::get_elements() = quadrature_dmn::get_elements();
    }

    {
      coarsegraining_domain<K_dmn, K>::get_size()     = quadrature_dmn::get_size();
      coarsegraining_domain<K_dmn, K>::get_weights()  = quadrature_dmn::get_weights();
      coarsegraining_domain<K_dmn, K>::get_elements() = quadrature_dmn::get_elements();
    }

    {
      coarsegraining_domain<K_dmn, K_PLUS_Q>::get_size()     = quadrature_dmn::get_size();
      coarsegraining_domain<K_dmn, K_PLUS_Q>::get_weights()  = quadrature_dmn::get_weights();
      coarsegraining_domain<K_dmn, K_PLUS_Q>::get_elements() = quadrature_dmn::get_elements();
    }

    {
      coarsegraining_domain<K_dmn, Q_MINUS_K>::get_size()     = quadrature_dmn::get_size();
      coarsegraining_domain<K_dmn, Q_MINUS_K>::get_weights()  = quadrature_dmn::get_weights();
      coarsegraining_domain<K_dmn, Q_MINUS_K>::get_elements() = quadrature_dmn::get_elements();
    }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::wannier_interpolation(int K_ind,
                                                                              FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& f_k,
                                                                              FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& f_q)
  {
    typedef interpolation_matrices<scalar_type, k_dmn_t, q_dmn_t> interpolation_matrices_type;

    if(not interpolation_matrices_type::is_initialized())
      interpolation_matrices_type::initialize(concurrency);

    LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& T = interpolation_matrices_type::get(K_ind);

    scalar_type alpha(1.);
    scalar_type beta (0.);

    scalar_type* A_ptr = &real(f_k(0));
    scalar_type* B_ptr = &T(0,0);
    scalar_type* C_ptr = &real(f_q(0));

    int M = 2*nu_nu::dmn_size();
    int K = k_dmn_t::dmn_size();
    int N = q_dmn_t::dmn_size();

    int LDA = 2*nu_nu::dmn_size();
    int LDB = T.get_global_size().first;
    int LDC = 2*nu_nu::dmn_size();

    LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'T', M, N, K, alpha, A_ptr, LDA, B_ptr, LDB, beta, C_ptr, LDC);
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename tmp_scalar_type, typename q_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::compute_I_q(tmp_scalar_type value,
                                                                    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q)
  {
    for(int q_ind=0; q_ind<q_dmn_t::dmn_size(); q_ind++)
      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          I_q(i,j,q_ind) = i==j? value : 0.;
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::compute_H_q(int K_ind,
                                                                    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& H_k,
                                                                    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q)
  {
    wannier_interpolation(K_ind, H_k, H_q);
  }

  /*****************************************
   ***                                   ***
   ***         Routines for DCA          ***
   ***                                   ***
   *****************************************/

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename q_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::compute_S_q(int K_ind, int w_ind,
                                                                    FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn  , w> >& S_K_w,
                                                                    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t   > >& S_q)
  {
    for(int q_ind=0; q_ind<q_dmn_t::dmn_size(); q_ind++)
      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          S_q(i,j,q_ind) = S_K_w(i,j,K_ind,w_ind);
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(int K_ind, int w_ind,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t   > >& H_k,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn  , w> >& S_K,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t   > >& I_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t   > >& H_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t   > >& S_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t   > >& G_q)
  {
    {
      std::complex<scalar_type> i_wm_min_mu;

      real(i_wm_min_mu) = parameters.get_chemical_potential();
      imag(i_wm_min_mu) = w::get_elements()[w_ind];

      compute_I_q(i_wm_min_mu, I_q);
    }

    compute_H_q(K_ind, H_k, H_q);

    compute_S_q(K_ind, w_ind, S_K, S_q);

    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> G_inv("G_inv", nu::dmn_size());

    for(int q_ind=0; q_ind<q_dmn_t::dmn_size(); q_ind++){

      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          G_inv(i,j) = I_q(i,j,q_ind)-H_q(i,j,q_ind)-S_q(i,j,q_ind);

      if(false)
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(G_inv);
      else
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(G_inv);

      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          G_q(i,j,q_ind) = G_inv(i,j);
    }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename q_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::compute_S_q(int K_ind, int w_ind,
                                                                    FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn  , w_REAL> >& S_K_w,
                                                                    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t        > >& S_q)
  {
    for(int q_ind=0; q_ind<q_dmn_t::dmn_size(); q_ind++)
      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          S_q(i,j,q_ind) = S_K_w(i,j,K_ind,w_ind);
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(int K_ind, int w_ind,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t        > >& H_k,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn  , w_REAL> >& S_K,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t        > >& I_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t        > >& H_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t        > >& S_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t        > >& G_q)
  {
    {
      std::complex<scalar_type> i_wm_min_mu;

      real(i_wm_min_mu) = w_REAL::get_elements()[w_ind]+parameters.get_chemical_potential();
      imag(i_wm_min_mu) = parameters.get_real_frequencies_off_set();

      compute_I_q(i_wm_min_mu, I_q);
    }

    compute_H_q(K_ind, H_k, H_q);

    compute_S_q(K_ind, w_ind, S_K, S_q);

    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> G_inv("G_inv", nu::dmn_size());

    for(int q_ind=0; q_ind<q_dmn_t::dmn_size(); q_ind++){

      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          G_inv(i,j) = I_q(i,j,q_ind)-H_q(i,j,q_ind)-S_q(i,j,q_ind);

      if(false)
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute(G_inv);
      else
        LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(G_inv);

      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          G_q(i,j,q_ind) = G_inv(i,j);
    }
  }

  /*
    template<typename parameters_type, typename K_dmn>
    template<typename scalar_type, typename k_dmn_t, typename w_dmn_t, typename tet_dmn_t>
    void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w_with_TIM(int K_ind, int w_ind,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t           > >& H_k,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn    , w_dmn_t> >& S_K,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t         > >& I_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t         > >& H_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t         > >& S_q,
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t         > >& G_q)
    {
    {
    std::complex<scalar_type> i_wm_min_mu;

    real(i_wm_min_mu) = parameters.get_chemical_potential();
    imag(i_wm_min_mu) = w_dmn_t::get_elements()[w_ind];

    compute_I_q(i_wm_min_mu, I_q);
    }

    compute_H_q(K_ind, H_k, H_q);

    compute_S_q(K_ind, w_ind, S_K, S_q);

    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> G_inv("G_inv", nu::dmn_size());

    for(int tet_ind=0; tet_ind<tet_dmn_t::dmn_size(); tet_ind++){

    for(int j=0; j<nu::dmn_size(); j++)
    for(int i=0; i<nu::dmn_size(); i++)
    G_inv(i,j) = I_q(i,j,tet_ind)-H_q(i,j,tet_ind)-S_q(i,j,tet_ind);

    if(false)
    LIN_ALG::GEINV<LIN_ALG::CPU>::execute(G_inv);
    else
    LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(G_inv);

    for(int j=0; j<nu::dmn_size(); j++)
    for(int i=0; i<nu::dmn_size(); i++)
    G_q(i,j,tet_ind) = G_inv(i,j);
    }
    }
  */

  /*
   *
   * p 122 AGD
   *
   */
  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_t(int K_ind, int t_ind,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& H_k,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q)
  {
    scalar_type f_val = 1;
    scalar_type t_val = t::get_elements()[t_ind];
    scalar_type beta  = parameters.get_beta();

    f_val = t_val<0? 1          : -1;
    t_val = t_val<0? t_val+beta : t_val;

    compute_I_q(parameters.get_chemical_potential(), I_q);

    compute_H_q(K_ind, H_k, H_q);

    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> H_m("H_m", nu::dmn_size());

    LIN_ALG::vector<scalar_type              , LIN_ALG::CPU> L("e_l", nu::dmn_size());
    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> V("V_l", nu::dmn_size());

    LIN_ALG::vector<scalar_type              , LIN_ALG::CPU> G_t("e_l", nu::dmn_size());

    G_q = 0.;

    for(int q_ind=0; q_ind<q_dmn_t::dmn_size(); q_ind++){

      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          H_m(i,j) = H_q(i,j,q_ind)-I_q(i,j,q_ind);

      if(false)
        LIN_ALG::GEEV<LIN_ALG::CPU>::execute('V', 'U', H_m, L, V);
      else
        LIN_ALG::GEEV<LIN_ALG::CPU>::execute_on_Greens_function_matrix('V', 'U', H_m, L, V);

      for(int i=0; i<nu::dmn_size(); i++){

        if(L[i]<0)
          G_t[i] = f_val*std::exp(L[i]*(beta-t_val))/(std::exp(L[i]*beta)+1.);
        else
          G_t[i] = f_val*std::exp(-L[i]*t_val)/(std::exp(-L[i]*beta)+1.);

        if(G_t[i]!=G_t[i]){
          cout << "\n\t warning in compute_G_q_t --> G_t[i] : " << G_t[i] << "\n";
          cout << "\n\tL[i] : " << L[i] << "\n";
          cout << "\n\tbeta : " << beta << "\n";
          cout << "\n\ttau  : " << t_val << "\n";
          cout << "\n\tstd::exp(L[i]*beta)  : "        << std::exp(L[i]*beta)         << "\n";
          cout << "\n\tstd::exp(L[i]*(beta-t_val)) : " << std::exp(L[i]*(beta-t_val)) << "\n";

          throw std::logic_error(__FUNCTION__);
        }

        //         G_t[i] = f_val*std::exp(L[i]*(beta-t_val))/(std::exp(L[i]*beta)+1.);
        //         if(G_t[i]!=G_t[i]){
        //           cout << "\n\t warning in compute_G_q_t --> G_t[i] : " << G_t[i] << "(L[i]*(beta-t_val) : " << L[i]*(beta-t_val) << ", L[i]*beta : " << L[i]*beta << ")\n";
        //         }
      }

      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          for(int l=0; l<nu::dmn_size(); l++)
            G_q(i,j,q_ind) += G_t[l]*real(conj(V(l,i))*V(l,j));
    }
  }

  /*****************************************
   ***                                   ***
   ***         Routines for DCA+         ***
   ***                                   ***
   *****************************************/

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::compute_S_q_from_A_k(int K_ind, int w_ind,
                                                                             FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& A_k,
                                                                             FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& A_q,
                                                                             FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& S_q)
  {
    wannier_interpolation(K_ind, A_k, A_q);

    scalar_type alpha = w::get_elements()[w_ind]>0 ? 1 : -1;

    transform_to_alpha::backward(alpha, S_q, A_q);
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename k_dmn_t, typename q_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(int K_ind, int w_ind,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& H_k,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn_t> >& A_k,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& I_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& H_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& A_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& S_q,
                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, q_dmn_t> >& G_q)
  {
    {
      std::complex<scalar_type> i_wm_min_mu;

      real(i_wm_min_mu) = parameters.get_chemical_potential();
      imag(i_wm_min_mu) = w::get_elements()[w_ind];

      compute_I_q(i_wm_min_mu, I_q);
    }

    compute_H_q(K_ind, H_k, H_q);

    compute_S_q_from_A_k(K_ind, w_ind, A_k, A_q, S_q);

    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> G_inv("G_inv", nu::dmn_size());

    for(int q_ind=0; q_ind<q_dmn_t::dmn_size(); q_ind++){

      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          G_inv(i,j) = I_q(i,j,q_ind)-H_q(i,j,q_ind)-S_q(i,j,q_ind);

      LIN_ALG::GEINV<LIN_ALG::CPU>::execute_on_Green_function_matrix(G_inv);

      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          G_q(i,j,q_ind) = G_inv(i,j);
    }
  }

  /***********************************************************************************
   ***                                                                             ***
   ***         tetrahedron-integration
   ***                                                                             ***
   ***********************************************************************************/

  /*
  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename tet_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::tetrahedron_integration_st(FUNC_LIB::function<             scalar_type , tet_dmn_t>&                 w_tet,
                                                                                   FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                                                                   FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int)
  {
    switch(DIMENSION)
      {
      case 1:
        tetrahedron_integration_st_1D(w_tet, G_tet, G_int);
        break;

      case 2:
        tetrahedron_integration_st_2D(w_tet, G_tet, G_int);
        break;

      case 3:
        tetrahedron_integration_st_3D(w_tet, G_tet, G_int);
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename tet_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::tetrahedron_integration_mt(FUNC_LIB::function<             scalar_type , tet_dmn_t>&                 w_tet,
                                                                                   FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                                                                   FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int)
  {
    switch(DIMENSION)
      {
      case 1:
        tetrahedron_integration_mt_1D(w_tet, G_tet, G_int);
        break;

      case 2:
        tetrahedron_integration_mt_2D(w_tet, G_tet, G_int);
        break;

      case 3:
        tetrahedron_integration_mt_3D(w_tet, G_tet, G_int);
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename tet_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::tetrahedron_integration_st_1D(FUNC_LIB::function<             scalar_type , tet_dmn_t>&                 w_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int)
  {
    int size = nu::dmn_size();

    for(int tet_ind=0; tet_ind<tet_dmn_t::dmn_size(); tet_ind += 2)
      {
        scalar_type volume = w_tet(tet_ind)+w_tet(tet_ind+1);

        std::complex<scalar_type>* G_0   = &G_tet(0,0,tet_ind+0);
        std::complex<scalar_type>* G_1   = &G_tet(0,0,tet_ind+1);
        std::complex<scalar_type>* G_ptr = &G_int(0,0);

        tetrahedron_integration_routines::execute(size, volume, G_0, G_1, G_ptr);
      }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename tet_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::tetrahedron_integration_mt_1D(FUNC_LIB::function<             scalar_type , tet_dmn_t>&                 w_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int)
  {
    int size = nu::dmn_size();

    for(int tet_ind=0; tet_ind<tet_dmn_t::dmn_size(); tet_ind += 2)
      {
        scalar_type volume = w_tet(tet_ind)+w_tet(tet_ind+1);

        std::complex<scalar_type>* G_0   = &G_tet(0,0,tet_ind+0);
        std::complex<scalar_type>* G_1   = &G_tet(0,0,tet_ind+1);
        std::complex<scalar_type>* G_ptr = &G_int(0,0);

        tetrahedron_integration_routines::execute(size, volume, G_0, G_1, G_ptr);
      }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename tet_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::tetrahedron_integration_st_2D(FUNC_LIB::function<             scalar_type , tet_dmn_t>&                 w_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int)
  {
    int size = nu::dmn_size();

    for(int tet_ind=0; tet_ind<tet_dmn_t::dmn_size(); )
      {
        scalar_type volume = w_tet(tet_ind)+w_tet(tet_ind+1)+w_tet(tet_ind+2);

        std::complex<scalar_type>* G_0   = &G_tet(0,0,tet_ind+0);
        std::complex<scalar_type>* G_1   = &G_tet(0,0,tet_ind+1);
        std::complex<scalar_type>* G_2   = &G_tet(0,0,tet_ind+2);
        std::complex<scalar_type>* G_ptr = &G_int(0,0);

        tetrahedron_integration_routines::execute(size, volume, G_0, G_1, G_2, G_ptr);

	tet_ind += 3;
      }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename tet_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::tetrahedron_integration_mt_2D(FUNC_LIB::function<scalar_type, tet_dmn_t>&                               w_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int)
  {
    int nr_threads = 8;

    tetrahedron_integration_functions<scalar_type> tetrahedron_integration_functions_obj;

    tetrahedron_integration_functions_obj.size = nu::dmn_size();

    tetrahedron_integration_functions_obj.w_tet_ptr = &w_tet;
    tetrahedron_integration_functions_obj.G_tet_ptr = &G_tet;

    tetrahedron_integration_functions_obj.G_int_vec.resize(nr_threads);

    COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY> parallelization_obj;

    parallelization_obj.execute(nr_threads, threaded_tetrahedron_integration_2D<scalar_type>, tetrahedron_integration_functions_obj);

    for(int j=0; j<nu::dmn_size(); j++)
      for(int i=0; i<nu::dmn_size(); i++)
        G_int(i,j) = 0;

    for(int l=0; l<nr_threads; l++)
      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          G_int(i,j) += tetrahedron_integration_functions_obj.G_int_vec[l](i,j);

    //     int size = nu::dmn_size();

    //     for(int tet_ind=0; tet_ind<tet_dmn_t::dmn_size(); tet_ind += 3)
    //       {
    //         scalar_type volume = w_tet(tet_ind)+w_tet(tet_ind+1)+w_tet(tet_ind+2);

    //         std::complex<scalar_type>* G_0   = &G_tet(0,0,tet_ind+0);
    //         std::complex<scalar_type>* G_1   = &G_tet(0,0,tet_ind+1);
    //         std::complex<scalar_type>* G_2   = &G_tet(0,0,tet_ind+2);
    //         std::complex<scalar_type>* G_ptr = &G_int(0,0);

    //         tetrahedron_integration_routines::execute(size, volume, G_0, G_1, G_2, G_ptr);
    //       }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type>
  void* coarsegraining_routines<parameters_type, K_dmn>::threaded_tetrahedron_integration_2D(void* void_ptr)
  {
    typedef tetrahedron_integration_functions<scalar_type> tetrahedron_functions_type;

    COMP_LIB::posix_data*       data_ptr      = static_cast<COMP_LIB::posix_data      *>(void_ptr);
    tetrahedron_functions_type* functions_ptr = static_cast<tetrahedron_functions_type*>(data_ptr->args);

    int id         = data_ptr->id;
    int nr_threads = data_ptr->nr_threads;

    int size = functions_ptr->size;

    FUNC_LIB::function<             scalar_type , tet_dmn_type>&                 w_tet = *(functions_ptr->w_tet_ptr);
    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_type> >& G_tet = *(functions_ptr->G_tet_ptr);
    FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&               G_int =  (functions_ptr->G_int_vec[id]);

    tet_dmn_type        tet_dmn;
    std::pair<int, int> tet_bounds = COMP_LIB::parallelization<COMP_LIB::POSIX_LIBRARY>::get_bounds(id, nr_threads, tet_dmn);

    for(int tet_ind=tet_bounds.first; tet_ind<tet_bounds.second; tet_ind++)
      {
        scalar_type volume = w_tet(tet_ind)+w_tet(tet_ind+1)+w_tet(tet_ind+2);

        std::complex<scalar_type>* G_0   = &G_tet(0,0,tet_ind+0);
        std::complex<scalar_type>* G_1   = &G_tet(0,0,tet_ind+1);
        std::complex<scalar_type>* G_2   = &G_tet(0,0,tet_ind+2);
        std::complex<scalar_type>* G_ptr = &G_int(0,0);

        tetrahedron_integration_routines::execute(size, volume, G_0, G_1, G_2, G_ptr);
      }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename tet_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::tetrahedron_integration_st_3D(FUNC_LIB::function<scalar_type, tet_dmn_t>&                 w_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int)
  {
    int size = nu::dmn_size();

    for(int tet_ind=0; tet_ind<tet_dmn_t::dmn_size(); tet_ind += 4)
      {
        scalar_type volume = w_tet(tet_ind)+w_tet(tet_ind+1)+w_tet(tet_ind+2)+w_tet(tet_ind+3);

        std::complex<scalar_type>* G_0   = &G_tet(0,0,tet_ind+0);
        std::complex<scalar_type>* G_1   = &G_tet(0,0,tet_ind+1);
        std::complex<scalar_type>* G_2   = &G_tet(0,0,tet_ind+2);
        std::complex<scalar_type>* G_3   = &G_tet(0,0,tet_ind+3);
        std::complex<scalar_type>* G_ptr = &G_int(0,0);

        tetrahedron_integration_routines::execute(size, volume, G_0, G_1, G_2, G_3, G_ptr);
      }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename scalar_type, typename tet_dmn_t>
  void coarsegraining_routines<parameters_type, K_dmn>::tetrahedron_integration_mt_3D(FUNC_LIB::function<scalar_type, tet_dmn_t>&                 w_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, tet_dmn_t> >& G_tet,
                                                                                      FUNC_LIB::function<std::complex<scalar_type>, dmn_2<nu, nu> >&            G_int)
  {
    int size = nu::dmn_size();

    for(int tet_ind=0; tet_ind<tet_dmn_t::dmn_size(); tet_ind += 4)
      {
        scalar_type volume = w_tet(tet_ind)+w_tet(tet_ind+1)+w_tet(tet_ind+2)+w_tet(tet_ind+3);

        std::complex<scalar_type>* G_0   = &G_tet(0,0,tet_ind+0);
        std::complex<scalar_type>* G_1   = &G_tet(0,0,tet_ind+1);
        std::complex<scalar_type>* G_2   = &G_tet(0,0,tet_ind+2);
        std::complex<scalar_type>* G_3   = &G_tet(0,0,tet_ind+3);
        std::complex<scalar_type>* G_ptr = &G_int(0,0);

        tetrahedron_integration_routines::execute(size, volume, G_0, G_1, G_2, G_3, G_ptr);
      }
  }
  */

}

#endif
