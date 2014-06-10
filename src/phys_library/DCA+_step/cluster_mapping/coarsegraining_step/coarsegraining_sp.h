//-*-C++-*-

#ifndef DCA_COARSEGRAINING_SP_H
#define DCA_COARSEGRAINING_SP_H

namespace DCA
{
  template<typename parameters_type, typename K_dmn>
  class coarsegraining_sp : public coarsegraining_routines<parameters_type, K_dmn>,
                            public tetrahedron_integration<parameters_type, K_dmn>
  {
#include "type_definitions.h"

    typedef typename K_dmn::parameter_type k_cluster_type;

    const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

    typedef typename parameters_type::concurrency_type concurrency_type;

#ifdef SINGLE_PRECISION_COARSEGRAINING
    typedef              float        scalar_type;
#else
    typedef              double       scalar_type;
#endif

    typedef std::complex<scalar_type> complex_type;

    typedef dmn_0<MATH_ALGORITHMS::tetrahedron_mesh<K_dmn> >             tetrahedron_dmn;

    typedef MATH_ALGORITHMS::gaussian_quadrature_domain<tetrahedron_dmn> quadrature_dmn;

    typedef dmn_0<coarsegraining_domain<K_dmn, K     > > q_dmn;
    typedef dmn_0<coarsegraining_domain<K_dmn, ORIGIN> > q_0_dmn;

    typedef dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_K     > > tet_dmn;
    typedef dmn_0<coarsegraining_domain<K_dmn, TETRAHEDRON_ORIGIN> > tet_0_dmn;

    typedef dmn_3<nu, nu, q_dmn  > nu_nu_q;
    typedef dmn_3<nu, nu, tet_dmn> nu_nu_tet;

  public:

    coarsegraining_sp(parameters_type& parameters_ref);

    ~coarsegraining_sp();

    void initialize();

    void reset_fine_q_mesh(int recursion, int rule, int period);

    template<typename other_scalar_type, typename k_dmn>
    void plot_H_q(function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn> >& H_0);

    template<typename other_scalar_type, typename k_dmn>
    void plot_S_q(function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w);

    template<typename other_scalar_type, typename r_dmn>
    void compute_phi_r(function<other_scalar_type, r_dmn>& phi_r);

    /*****************************************
     ***                                   ***
     ***         Routines for DCA          ***
     ***                                   ***
     *****************************************/

    template<typename other_scalar_type>
    void compute_S_K_w(function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& S_k_w,
                       function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& S_K_w);

    template<typename other_scalar_type, typename k_dmn>
    void compute_G_K_w(function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn> >&    H_0,
                       function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& S_K_w,
                       function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& G_K_w);

    template<typename other_scalar_type, typename k_dmn, typename w_dmn>
    void compute_G_K_w_with_TIM(function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn> >&        H_0,
                                function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn> >& S_K_w,
                                function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn> >& G_K_w);

    /*****************************************
     ***                                   ***
     ***         Routines for DCA+         ***
     ***                                   ***
     *****************************************/

    template<typename other_scalar_type, typename k_dmn>
    void compute_S_K_w(function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w,
                       function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& S_K_w);


    template<typename other_scalar_type, typename k_dmn>
    void compute_G_K_w(function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn> >&    H_0,
                       function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w,
                       function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& G_K_w);

    template<typename other_scalar_type, typename k_dmn, typename w_dmn>
    void compute_G_K_w_with_TIM(function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn> >&        H_0,
                                function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w_dmn> >& S_K_w,
                                function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn> >& G_K_w);


    template<typename other_scalar_type, typename k_dmn>
    void compute_G0_K_t(function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn   > >& H_0,
                        function<             other_scalar_type , dmn_4<nu, nu, K_dmn, t> >& G0_k_w);

  private:

    void update_shell(int i, int N);

    template<typename other_scalar_type, typename r_dmn>
    void plot_phi_r(function<other_scalar_type, r_dmn>& phi_r);

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;

    bool              initialized;

    // tetrahedron q-points
    function<             scalar_type,  tet_dmn  > w_tet;

    function<std::complex<scalar_type>, nu_nu_tet> I_tet;
    function<std::complex<scalar_type>, nu_nu_tet> H_tet;
    function<std::complex<scalar_type>, nu_nu_tet> S_tet;
    function<std::complex<scalar_type>, nu_nu_tet> A_tet;
    function<std::complex<scalar_type>, nu_nu_tet> G_tet;

    // gaussian q-points
    function<             scalar_type,  q_dmn  > w_q;

    function<std::complex<scalar_type>, nu_nu_q> I_q;
    function<std::complex<scalar_type>, nu_nu_q> H_q;
    function<std::complex<scalar_type>, nu_nu_q> S_q;
    function<std::complex<scalar_type>, nu_nu_q> A_q;
    function<std::complex<scalar_type>, nu_nu_q> G_q;
  };

  template<typename parameters_type, typename K_dmn>
  coarsegraining_sp<parameters_type, K_dmn>::coarsegraining_sp(parameters_type& parameters_ref):
    coarsegraining_routines<parameters_type, K_dmn>(parameters_ref),
    tetrahedron_integration<parameters_type, K_dmn>(parameters_ref),

    parameters (parameters_ref),
    concurrency(parameters.get_concurrency()),

    initialized(false),

    // tetrahedron q-points
    w_tet("w_tet"),

    I_tet("I_tet"),
    H_tet("H_tet"),
    S_tet("S_tet"),
    A_tet("A_tet"),
    G_tet("G_tet"),

    // gaussian q-points
    w_q("w_q"),

    I_q("I_q"),
    H_q("H_q"),
    S_q("S_q"),
    A_q("A_q"),
    G_q("G_q")
  {
    initialize();
  }

  template<typename parameters_type, typename K_dmn>
  coarsegraining_sp<parameters_type, K_dmn>::~coarsegraining_sp()
  {}

  template<typename parameters_type, typename K_dmn>
  void coarsegraining_sp<parameters_type, K_dmn>::update_shell(int i, int N)
  {
    int tmp = i;

    if( concurrency.id() == concurrency.first() && N > 10 && (tmp % (N/10)) == 0 )
      {
        std::stringstream ss;

        ss << scientific;
        ss.precision(6);

        ss << "\t\t\t"<< double(i)/double(N)*100. << " % completed \t ";
        ss << print_time();
	ss << "\n";

        cout << ss.str();
      }
  }

  template<typename parameters_type, typename K_dmn>
  void coarsegraining_sp<parameters_type, K_dmn>::initialize()
  {
    //     if(concurrency.id() == concurrency.first())
    //       cout << "\n\n\t initialization of coarsegraining_sp has started | time = " << print_time();

    //SHOW::plot_points(K_dmn::get_elements());

    {
      this->compute_tetrahedron_mesh(parameters.get_k_mesh_refinement(),
                                     parameters.get_number_of_periods());

      this->compute_gaussian_mesh(parameters.get_k_mesh_refinement(),
                                  parameters.get_gaussian_quadrature_rule(),
                                  parameters.get_number_of_periods());
    }

    {
      w_q  .reset();
      w_tet.reset();

      //       for(int l=0; l<w_q.size(); l++)
      //         w_q(l) = quadrature_dmn::get_weights()[l];

      for(int l=0; l<w_q.size(); l++)
        w_q(l) = q_dmn::parameter_type::get_weights()[l];

      for(int l=0; l<w_tet.size(); l++)
        w_tet(l) = tet_dmn::parameter_type::get_weights()[l];
    }

    {
      I_tet.reset();
      H_tet.reset();
      S_tet.reset();
      A_tet.reset();
      G_tet.reset();

      I_q.reset();
      H_q.reset();
      S_q.reset();
      A_q.reset();
      G_q.reset();
    }

    interpolation_matrices<scalar_type, k_HOST, q_dmn>::initialize(concurrency);

    //     if(concurrency.id() == concurrency.first())
    //       cout << "\t initialization of coarsegraining_sp has ended   | time = " << print_time() << "\n";
  }

  template<typename parameters_type, typename K_dmn>
  template<typename other_scalar_type, typename k_dmn>
  void coarsegraining_sp<parameters_type, K_dmn>::plot_H_q(function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn> >& H_0)
  {
    /*
      function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn> > H_k("H_k");

      for(int k_ind=0; k_ind<k_dmn::dmn_size(); k_ind++)
      for(int j=0; j<nu::dmn_size(); j++)
      for(int i=0; i<nu::dmn_size(); i++)
      H_k(i,j,k_ind) = H_0(i,j,k_ind);

      std::vector<scalar_type> x(q_dmn::dmn_size()*K_dmn::dmn_size());
      std::vector<scalar_type> y(q_dmn::dmn_size()*K_dmn::dmn_size());
      std::vector<scalar_type> z(q_dmn::dmn_size()*K_dmn::dmn_size());

      for(int K_ind=0; K_ind<K_dmn::dmn_size(); K_ind++){

      compute_H_q(K_ind, H_k);

      for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++){
      x[q_ind+K_ind*q_dmn::dmn_size()] = q_0_dmn::parameter_type::get_elements()[q_ind][0]+K_dmn::get_elements()[K_ind][0];
      y[q_ind+K_ind*q_dmn::dmn_size()] = q_0_dmn::parameter_type::get_elements()[q_ind][1]+K_dmn::get_elements()[K_ind][1];
      z[q_ind+K_ind*q_dmn::dmn_size()] = real(H_q(0,0,q_ind));
      }
      }

      SHOW::plot_points(x,y);
      SHOW::heatmap(x,y,z);
    */
  }

  template<typename parameters_type, typename K_dmn>
  template<typename other_scalar_type, typename k_dmn>
  void coarsegraining_sp<parameters_type, K_dmn>::plot_S_q(function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w)
  {
    /*
      std::vector<scalar_type> x(q_dmn::dmn_size()*K_dmn::dmn_size());
      std::vector<scalar_type> y(q_dmn::dmn_size()*K_dmn::dmn_size());
      std::vector<scalar_type> z(q_dmn::dmn_size()*K_dmn::dmn_size());

      function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn> > S_k;
      function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn> > A_k;

      for(int k_ind=0; k_ind<k_dmn::dmn_size(); k_ind++)
      for(int j=0; j<nu::dmn_size(); j++)
      for(int i=0; i<nu::dmn_size(); i++)
      S_k(i,j,k_ind) = S_k_w(i,j,k_ind,w::dmn_size()/2);

      transform_to_alpha::forward(scalar_type(1), S_k, A_k);

      for(int K_ind=0; K_ind<K_dmn::dmn_size(); K_ind++){

      compute_S_q_from_A_k(K_ind, w::dmn_size()/2, A_k);

      for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++){
      x[q_ind+K_ind*q_dmn::dmn_size()] = q_0_dmn::parameter_type::get_elements()[q_ind][0]+K_dmn::get_elements()[K_ind][0];
      y[q_ind+K_ind*q_dmn::dmn_size()] = q_0_dmn::parameter_type::get_elements()[q_ind][1]+K_dmn::get_elements()[K_ind][1];
      z[q_ind+K_ind*q_dmn::dmn_size()] = imag(S_q(0,0,q_ind));
      }
      }

      SHOW::plot_points(x,y);
      SHOW::heatmap(x,y,z);
    */
  }

  template<typename parameters_type, typename K_dmn>
  template<typename other_scalar_type, typename r_dmn>
  void coarsegraining_sp<parameters_type, K_dmn>::plot_phi_r(function<other_scalar_type, r_dmn>& phi_r)
  {
    std::vector<double> x;
    std::vector<double> y;

    std::vector<std::vector<double> > super_basis = r_dmn::parameter_type::get_super_basis_vectors();
    for(int r_ind=0; r_ind<r_dmn::dmn_size(); r_ind++){

      std::vector<double>               r_vec  = r_dmn::get_elements()[r_ind];
      std::vector<std::vector<double> > r_vecs = cluster_operations::equivalent_vectors(r_vec, super_basis);

      x.push_back(std::sqrt(r_vecs[0][0]*r_vecs[0][0]+r_vecs[0][1]*r_vecs[0][1]));
      y.push_back((phi_r(r_ind)));
    }

    SHOW::plot_points(x, y);
  }

  template<typename parameters_type, typename K_dmn>
  template<typename other_scalar_type, typename r_dmn>
  void coarsegraining_sp<parameters_type, K_dmn>::compute_phi_r(function<other_scalar_type, r_dmn>& phi_r)
  {
    MATH_ALGORITHMS::tetrahedron_mesh<k_cluster_type> mesh(parameters.get_k_mesh_refinement());

    quadrature_dmn::translate_according_to_period(parameters.get_number_of_periods(), mesh);

    std::vector<MATH_ALGORITHMS::tetrahedron<DIMENSION> >& tetrahedra = mesh.get_tetrahedra();

    {
      phi_r = 0.;

      r_dmn r_domain;
      std::pair<int, int> bounds = concurrency.get_bounds(r_domain);

      std::vector<std::vector<double> > super_basis = r_dmn::parameter_type::get_super_basis_vectors();

      for(int l=bounds.first; l<bounds.second; l++)
        {
          std::vector<double>               r_vec  = r_dmn::get_elements()[l];
          std::vector<std::vector<double> > r_vecs = cluster_operations::equivalent_vectors(r_vec, super_basis);

          //      for(int tet_ind=0; tet_ind<tetrahedra.size(); tet_ind++)
          //        {
          //          std::complex<double> val_0 = tetrahedra[tet_ind].integrate_harmonic(r_vecs[0]);
          //          std::complex<double> val_1 = tetrahedron_routines_harmonic_function::execute(r_vecs[0],
          //                                                                                       tetrahedra[tet_ind]);

          //          if(abs(val_0-val_1)>1.e-6){
          //            VECTOR_OPERATIONS::PRINT(r_vecs[0]);
          //            cout << "\t" << __FUNCTION__ << "\t" << val_0 << "\t" << val_1 << " be carefull !! \n";
          //          }
          //          else
          //            cout << "\t" << __FUNCTION__ << "\t" << val_0 << "\t" << val_1 << " OK !! \n";
          //        }

          //      if(false)
          //        {
          //          for(int r_ind=0; r_ind<r_vecs.size(); r_ind++)
          //            for(int tet_ind=0; tet_ind<tetrahedra.size(); tet_ind++)
          //              phi_r(l) += real(tetrahedra[tet_ind].integrate_harmonic(r_vecs[r_ind]))/r_vecs.size();
          //        }
          //      else
          //        {
          for(int r_ind=0; r_ind<r_vecs.size(); r_ind++)
            for(int tet_ind=0; tet_ind<tetrahedra.size(); tet_ind++)
              phi_r(l) += real(tetrahedron_routines_harmonic_function::execute(r_vecs[0],
                                                                               tetrahedra[tet_ind]))/r_vecs.size();
          //        }
        }

      concurrency.sum(phi_r);

      {
        scalar_type V_K=0;
        for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
          V_K += w_q(q_ind);

        phi_r /= V_K;
      }
    }
  }

  /*****************************************
   ***                                   ***
   ***         Routines for DCA          ***
   ***                                   ***
   *****************************************/

  template<typename parameters_type, typename K_dmn>
  template<typename other_scalar_type>
  void coarsegraining_sp<parameters_type, K_dmn>::compute_S_K_w(function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& S_k_w,
                                                                function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& S_K_w)
  {
    for(int l=0; l<S_k_w.size(); l++)
      S_K_w(l) = S_k_w(l);
  }

  template<typename parameters_type, typename K_dmn>
  template<typename other_scalar_type, typename k_dmn>
  void coarsegraining_sp<parameters_type, K_dmn>::compute_G_K_w(function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn> >&    H_0,
                                                                function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& S_K_w,
                                                                function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& G_K_w)
  {
    //cout << "\n\n\n\t" << __FUNCTION__ << "\n\n\n";

    function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn> >    H_k("H_k");
    function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn, w> > S_K("S_K");

    for(int l=0; l<S_K_w.size(); l++)
      S_K(l) = S_K_w(l);

    for(int k_ind=0; k_ind<k_dmn::dmn_size(); k_ind++)
      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          H_k(i,j,k_ind) = H_0(i,j,k_ind);

    G_K_w = 0.;

    dmn_2<K_dmn, w> K_wm_dmn;
    std::pair<int, int> bounds = concurrency.get_bounds(K_wm_dmn);

    int coor[2];
    for(int l=bounds.first; l<bounds.second; l++)
      {
        K_wm_dmn.linind_2_subind(l, coor);

        this->compute_G_q_w(coor[0], coor[1], H_k, S_K, I_q, H_q, S_q, G_q);

        for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
          for(int j=0; j<nu::dmn_size(); j++)
            for(int i=0; i<nu::dmn_size(); i++)
              G_K_w(i,j,coor[0],coor[1]) += G_q(i,j,q_ind)*w_q(q_ind);
      }

    concurrency.sum(G_K_w);

    {
      scalar_type V_K=0;
      for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
        V_K += w_q(q_ind);

      G_K_w /= V_K;
    }

    //SHOW::execute_on_bands(G_K_w);
    //throw std::logic_error(__FUNCTION__);
  }

  template<typename parameters_type, typename K_dmn>
  template<typename other_scalar_type, typename k_dmn, typename w_dmn>
  void coarsegraining_sp<parameters_type, K_dmn>::compute_G_K_w_with_TIM(function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn> >&        H_0,
                                                                         function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn> >& S_K_w,
                                                                         function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn> >& G_K_w)
  {
//     if(concurrency.id()==0)
//       cout << "\n\t start " << __FUNCTION__ << " : " << print_time() << endl;

    function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn       > > H_k("H_k");
    function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn, w_dmn> > S_K("S_K");

    for(int l=0; l<S_K_w.size(); l++)
      S_K(l) = S_K_w(l);

    for(int k_ind=0; k_ind<k_dmn::dmn_size(); k_ind++)
      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          H_k(i,j,k_ind) = H_0(i,j,k_ind);

    G_K_w = 0.;

    /*
    {
      for(int tet_ind=0; tet_ind<tet_dmn::dmn_size(); tet_ind += 4)
	{
	  VECTOR_OPERATIONS::PRINT(tet_dmn::get_elements()[tet_ind+0]); cout << "\n";
	  VECTOR_OPERATIONS::PRINT(tet_dmn::get_elements()[tet_ind+1]); cout << "\n";
	  VECTOR_OPERATIONS::PRINT(tet_dmn::get_elements()[tet_ind+2]); cout << "\n";
	  VECTOR_OPERATIONS::PRINT(tet_dmn::get_elements()[tet_ind+3]); cout << "\n";
	}
      cout << "\n";

      for(int tet_ind=0; tet_ind<tet_dmn::dmn_size(); tet_ind += 4)
	cout << tet_ind+0 << "\t" << tet_ind+1 << "\t" << tet_ind+2 << "\t" << tet_ind+3 << "\n";
      cout << "\n";

      assert(false);
    }
    */

    dmn_2<K_dmn, w_dmn> K_wm_dmn;
    std::pair<int, int> bounds = concurrency.get_bounds(K_wm_dmn);

    int coor[2];
    for(int l=bounds.first; l<bounds.second; l++)
      {
//         update_shell(l-bounds.first, bounds.second-bounds.first+1);

        K_wm_dmn.linind_2_subind(l, coor);

        this->compute_G_q_w(coor[0], coor[1], H_k, S_K, I_tet, H_tet, S_tet, G_tet);

//         if(l==w_dmn::dmn_size()/2-10 or l==w_dmn::dmn_size()/2+10)
//           {
//             //cout << w_dmn::get_elements()[l] << "\t";
//           }

        {
          function<std::complex<scalar_type>, dmn_2<nu, nu> > G_int;

          //this->tetrahedron_integration_st(w_tet, G_tet, G_int);
          this->tetrahedron_integration_mt(w_tet, G_tet, G_int);

          for(int j=0; j<nu::dmn_size(); j++)
            for(int i=0; i<nu::dmn_size(); i++)
              G_K_w(i,j,coor[0],coor[1]) += G_int(i,j);
        }
      }

    concurrency.sum(G_K_w);

    {
      scalar_type V_K=0;
      for(int tet_ind=0; tet_ind<tet_dmn::dmn_size(); tet_ind++)
        V_K += w_tet(tet_ind);

      G_K_w /= V_K;
    }

//     if(concurrency.id()==0)
//       cout << "\n\t end   " << __FUNCTION__ << " : " << print_time() << endl;
  }

  /*****************************************
   ***                                   ***
   ***         Routines for DCA+         ***
   ***                                   ***
   *****************************************/

  template<typename parameters_type, typename K_dmn>
  template<typename other_scalar_type, typename k_dmn>
  void coarsegraining_sp<parameters_type, K_dmn>::compute_S_K_w(function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w,
                                                                function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& S_K_w)
  {
    S_K_w = 0.;

    function<std::complex<      scalar_type>, dmn_3<nu, nu, k_dmn   > > A_k("A_k");
    function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w> > A_k_w("A_k_w");

    transform_to_alpha::forward (1., S_k_w, A_k_w);

    dmn_2<K_dmn, w> K_wm_dmn;
    std::pair<int, int> bounds = concurrency.get_bounds(K_wm_dmn);

    int coor[2];
    for(int l=bounds.first; l<bounds.second; l++)
      {
        K_wm_dmn.linind_2_subind(l, coor);

        for(int k_ind=0; k_ind<k_dmn::dmn_size(); k_ind++)
          for(int j=0; j<nu::dmn_size(); j++)
            for(int i=0; i<nu::dmn_size(); i++)
              A_k(i,j,k_ind) = A_k_w(i,j,k_ind,coor[1]);

        this->compute_S_q_from_A_k(coor[0], coor[1], A_k, A_q, S_q);

        for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
          for(int j=0; j<nu::dmn_size(); j++)
            for(int i=0; i<nu::dmn_size(); i++)
              S_K_w(i,j,coor[0],coor[1]) += S_q(i,j,q_ind)*w_q(q_ind);
      }

    concurrency.sum(S_K_w);

    {
      scalar_type V_K=0;
      for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
        V_K += w_q(q_ind);

      S_K_w /= V_K;
    }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename other_scalar_type, typename k_dmn>
  void coarsegraining_sp<parameters_type, K_dmn>::compute_G_K_w(function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn> >&    H_0,
                                                                function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w,
                                                                function<std::complex<other_scalar_type>, dmn_4<nu, nu, K_dmn, w> >& G_K_w)
  {
    G_K_w = 0.;

    function<std::complex<      scalar_type>, dmn_3<nu, nu, k_dmn> > H_k("H_k");
    function<std::complex<      scalar_type>, dmn_3<nu, nu, k_dmn> > A_k("A_k");

    function<std::complex<other_scalar_type>, dmn_4<nu, nu, k_dmn, w> > A_k_w("A_k_w");

    transform_to_alpha::forward(1., S_k_w, A_k_w);

    for(int k_ind=0; k_ind<k_dmn::dmn_size(); k_ind++)
      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          H_k(i,j,k_ind) = H_0(i,j,k_ind);

    dmn_2<K_dmn, w> K_wm_dmn;
    std::pair<int, int> bounds = concurrency.get_bounds(K_wm_dmn);

    int coor[2];
    for(int l=bounds.first; l<bounds.second; l++)
      {
        K_wm_dmn.linind_2_subind(l, coor);

        for(int k_ind=0; k_ind<k_dmn::dmn_size(); k_ind++)
          for(int j=0; j<nu::dmn_size(); j++)
            for(int i=0; i<nu::dmn_size(); i++)
              A_k(i,j,k_ind) = A_k_w(i,j,k_ind,coor[1]);

        this->compute_G_q_w(coor[0], coor[1], H_k, A_k, I_q, H_q, A_q, S_q, G_q);

        for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
          for(int j=0; j<nu::dmn_size(); j++)
            for(int i=0; i<nu::dmn_size(); i++)
              G_K_w(i,j,coor[0],coor[1]) += G_q(i,j,q_ind)*w_q(q_ind);
      }

    concurrency.sum(G_K_w);

    {
      scalar_type V_K=0;
      for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
        V_K += w_q(q_ind);

      G_K_w /= V_K;
    }
  }

  /*
    template<typename parameters_type, typename K_dmn>
    template<typename k_dmn>
    void coarsegraining_sp<parameters_type, K_dmn>::compute_G_q_t(int K_ind, int t_ind,
    function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn> >& H_k)
    {
    scalar_type f_val = 1;
    scalar_type t_val = t::get_elements()[t_ind];
    scalar_type beta  = parameters.get_beta();

    f_val = t_val<0? 1          : -1;
    t_val = t_val<0? t_val+beta : t_val;

    //compute_I_q(parameters.get_chemical_potential());
    this->compute_I_q(parameters.get_chemical_potential(), I_q);

    //compute_H_q(K_ind, H_0);
    this->compute_H_q(K_ind, H_k, H_q);

    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> H_m("H_m", nu::dmn_size());

    LIN_ALG::vector<scalar_type              , LIN_ALG::CPU> L("e_l", nu::dmn_size());
    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> V("V_l", nu::dmn_size());

    LIN_ALG::vector<scalar_type              , LIN_ALG::CPU> G_t("e_l", nu::dmn_size());

    G_q = 0.;

    for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++){

    for(int j=0; j<nu::dmn_size(); j++)
    for(int i=0; i<nu::dmn_size(); i++)
    H_m(i,j) = H_q(i,j,q_ind)-I_q(i,j,q_ind);

    LIN_ALG::GEEV<LIN_ALG::CPU>::execute('V', 'U', H_m, L, V);

    for(int i=0; i<nu::dmn_size(); i++)
    G_t[i] = f_val*std::exp(L[i]*(beta-t_val))/(std::exp(L[i]*beta)+1.);

    for(int j=0; j<nu::dmn_size(); j++)
    for(int i=0; i<nu::dmn_size(); i++)
    for(int l=0; l<nu::dmn_size(); l++)
    G_q(i,j,q_ind) += G_t[l]*real(conj(V(l,i))*V(l,j));
    }
    }
  */

  template<typename parameters_type, typename K_dmn>
  template<typename other_scalar_type, typename k_dmn>
  void coarsegraining_sp<parameters_type, K_dmn>::compute_G0_K_t(function<std::complex<other_scalar_type>, dmn_3<nu, nu, k_dmn   > >& H_0,
                                                                 function<             other_scalar_type , dmn_4<nu, nu, K_dmn, t> >& G_K_t)
  {
    function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn> > H_k("H_k");

    for(int k_ind=0; k_ind<k_dmn::dmn_size(); k_ind++)
      for(int j=0; j<nu::dmn_size(); j++)
        for(int i=0; i<nu::dmn_size(); i++)
          H_k(i,j,k_ind) = H_0(i,j,k_ind);

    G_K_t = 0.;

    dmn_2<K_dmn, t> K_t_dmn;
    std::pair<int, int> bounds = concurrency.get_bounds(K_t_dmn);

    int coor[2];
    for(int l=bounds.first; l<bounds.second; l++)
      {
        K_t_dmn.linind_2_subind(l, coor);

        this->compute_G_q_t(coor[0], coor[1], H_k, I_q, H_q, G_q);

        for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
          for(int j=0; j<nu::dmn_size(); j++)
            for(int i=0; i<nu::dmn_size(); i++)
              G_K_t(i,j,coor[0],coor[1]) += real(G_q(i,j,q_ind))*w_q(q_ind);
      }

    concurrency.sum(G_K_t);

    {
      scalar_type V_K=0;
      for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
        V_K += w_q(q_ind);

      G_K_t /= V_K;
    }
  }

}

#endif
