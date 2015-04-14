//-*-C++-*-

#ifndef DCA_COARSEGRAINING_TP_H
#define DCA_COARSEGRAINING_TP_H

namespace DCA
{

  template<typename parameters_type, typename K_dmn>
  class coarsegraining_tp : public coarsegraining_routines<parameters_type, K_dmn>
  {
#include "type_definitions.h"

    typedef typename K_dmn::parameter_type k_cluster_type;

    const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

    typedef typename parameters_type::profiler_type    profiler_type;
    typedef typename parameters_type::concurrency_type concurrency_type;

    //typedef              float        scalar_type;
    typedef              double       scalar_type;

    typedef std::complex   <scalar_type>               complex_type;
    typedef LIN_ALG::matrix<scalar_type, LIN_ALG::CPU> matrix_type;

    typedef dmn_0<MATH_ALGORITHMS::tetrahedron_mesh<K_dmn> >             tetrahedron_dmn;
    typedef MATH_ALGORITHMS::gaussian_quadrature_domain<tetrahedron_dmn> quadrature_dmn;

    typedef dmn_0<coarsegraining_domain<K_dmn, K        > > q_dmn;
    typedef dmn_0<coarsegraining_domain<K_dmn, K_PLUS_Q > > q_plus_Q_dmn;
    typedef dmn_0<coarsegraining_domain<K_dmn, Q_MINUS_K> > Q_min_q_dmn;

    typedef dmn_3<nu, nu, q_dmn>        nu_nu_q;
    typedef dmn_3<nu, nu, q_plus_Q_dmn> nu_nu_q_plus_Q;
    typedef dmn_3<nu, nu, Q_min_q_dmn>  nu_nu_Q_min_q;

  public:

    coarsegraining_tp(parameters_type& parameters_ref);

    ~coarsegraining_tp();

    void initialize();

    template<typename w_dmn_t>
    void execute(FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu , nu , k_HOST         > >& H_k,
                 FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu , nu , k_HOST, w      > >& Sigma,
                 FUNC_LIB::function<std::complex<scalar_type>, dmn_4<b_b, b_b, K_dmn, w_dmn_t> >& chi);

    template<typename w_dmn_t>
    void plot(FUNC_LIB::function<std::complex<scalar_type>, dmn_4<b_b, b_b, K_dmn, w_dmn_t> >& chi);
    
  private:

    template<typename w_dmn_t>
    void compute_tp(FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu , nu , k_HOST         > >& H_k,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu , nu , k_HOST, w      > >& Sigma,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_4<b_b, b_b, K_dmn , w_dmn_t> >& chi);

    template<typename w_dmn_t>
    void compute_phi(FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu , nu , k_HOST         > >& H_k,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu , nu , k_HOST, w      > >& S_k_w,
                     FUNC_LIB::function<std::complex<scalar_type>, dmn_4<b_b, b_b, K_dmn , w_dmn_t> >& phi);


    void find_w1_and_w2(std::vector<double>& elements, int& w_ind, int& w1, int& w2);

    void compute_bubble(FUNC_LIB::function<std::complex<scalar_type>, dmn_3<b_b, b_b, q_dmn> >& bubble);

    double get_integration_factor();

  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;

    FUNC_LIB::function<             scalar_type,  q_dmn  > w_q;

    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q> I_q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q> H_q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q> S_q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q> A_q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q> G_q;

    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q_plus_Q> I_q_plus_Q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q_plus_Q> H_q_plus_Q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q_plus_Q> S_q_plus_Q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q_plus_Q> A_q_plus_Q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_q_plus_Q> G_q_plus_Q;

    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_Q_min_q> I_Q_min_q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_Q_min_q> H_Q_min_q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_Q_min_q> S_Q_min_q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_Q_min_q> A_Q_min_q;
    FUNC_LIB::function<std::complex<scalar_type>, nu_nu_Q_min_q> G_Q_min_q;

    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<b_b, b_b, q_dmn> > bubble_q;
  };

  template<typename parameters_type, typename K_dmn>
  coarsegraining_tp<parameters_type, K_dmn>::coarsegraining_tp(parameters_type& parameters_ref):
    coarsegraining_routines<parameters_type, K_dmn>(parameters_ref),

    parameters (parameters_ref),
    concurrency(parameters.get_concurrency()),

    w_q("w_q"),

    I_q("I_q"),
    H_q("H_q"),
    S_q("S_q"),
    A_q("A_q"),
    G_q("G_q"),

    I_q_plus_Q("I_q_plus_Q"),
    H_q_plus_Q("H_q_plus_Q"),
    S_q_plus_Q("S_q_plus_Q"),
    A_q_plus_Q("A_q_plus_Q"),
    G_q_plus_Q("G_q_plus_Q"),

    I_Q_min_q("I_Q_min_q"),
    H_Q_min_q("H_Q_min_q"),
    S_Q_min_q("S_Q_min_q"),
    A_Q_min_q("A_Q_min_q"),
    G_Q_min_q("G_Q_min_q"),

    bubble_q("bubble_q")
  {
    initialize();
  }

  template<typename parameters_type, typename K_dmn>
  coarsegraining_tp<parameters_type, K_dmn>::~coarsegraining_tp()
  {}

  template<typename parameters_type, typename K_dmn>
  void coarsegraining_tp<parameters_type, K_dmn>::initialize()
  {
    //SHOW::plot_points(K_dmn::get_elements());
    
    {
      if(true)
	coarsegraining_routines<parameters_type, K_dmn>::compute_gaussian_mesh(parameters.get_k_mesh_refinement(),
									       parameters.get_gaussian_quadrature_rule(),
									       parameters.get_number_of_periods());
      else
	coarsegraining_routines<parameters_type, K_dmn>::compute_gaussian_mesh(1, 0, 0);
    }

    {
      w_q.reset();

      for(int l=0; l<w_q.size(); l++)
        w_q(l) = quadrature_dmn::get_weights()[l];
    }

    {
      I_q.reset();
      H_q.reset();
      S_q.reset();
      A_q.reset();
      G_q.reset();
    }

    {
      I_q_plus_Q.reset();
      H_q_plus_Q.reset();
      S_q_plus_Q.reset();
      A_q_plus_Q.reset();
      G_q_plus_Q.reset();
    }

    {
      I_Q_min_q.reset();
      H_Q_min_q.reset();
      S_Q_min_q.reset();
      A_Q_min_q.reset();
      G_Q_min_q.reset();
    }

    bubble_q.reset();
  }

  template<typename parameters_type, typename K_dmn>
  template<typename w_dmn_t>
  void coarsegraining_tp<parameters_type, K_dmn>::execute(FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu , nu , k_HOST         > >& H_k,
                                                           FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu , nu , k_HOST, w      > >& Sigma,
                                                           FUNC_LIB::function<std::complex<scalar_type>, dmn_4<b_b, b_b, K_dmn , w_dmn_t> >& chi)
  {
    int Q_ind = cluster_operations::index(parameters.get_q_vector(), K_dmn::get_elements());

    switch(parameters.get_vertex_measurement_type())
      {
      case PARTICLE_HOLE_CHARGE:
      case PARTICLE_HOLE_MAGNETIC:
      case PARTICLE_HOLE_TRANSVERSE:
        {
          interpolation_matrices<scalar_type, k_HOST, q_dmn       >::initialize(concurrency);
          interpolation_matrices<scalar_type, k_HOST, q_plus_Q_dmn>::initialize(concurrency, Q_ind);

          compute_tp(H_k, Sigma, chi);
        }
        break;

      case PARTICLE_PARTICLE_SUPERCONDUCTING:
        {
          interpolation_matrices<scalar_type, k_HOST, q_dmn      >::initialize(concurrency);
          interpolation_matrices<scalar_type, k_HOST, Q_min_q_dmn>::initialize(concurrency, Q_ind);

          compute_phi(H_k, Sigma, chi);
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }
  
  template<typename parameters_type, typename K_dmn>
  template<typename w_dmn_t>
  void coarsegraining_tp<parameters_type, K_dmn>::plot(FUNC_LIB::function<std::complex<scalar_type>, dmn_4<b_b, b_b, K_dmn, w_dmn_t> >& chi)
  {
    {
      Gnuplot plot_obj;

      FUNC_LIB::function<scalar_type, w_dmn_t> phi_w("phi_w");

      for(int m2=0; m2<b::dmn_size(); m2++){
	for(int m1=0; m1<b::dmn_size(); m1++){
	  for(int n2=0; n2<b::dmn_size(); n2++){
	    for(int n1=0; n1<b::dmn_size(); n1++){
	      for(int k_ind=0; k_ind<K_dmn::dmn_size(); k_ind++){
		
		if(abs(K_dmn::get_elements()[k_ind][0]+K_dmn::get_elements()[k_ind][1]-M_PI)<1.e-6){
		  for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++)
		    phi_w(w_ind) = real(chi(n1,n2,m1,m2,k_ind,w_ind));

		  SHOW::execute(plot_obj, phi_w);
		}
	      }
	    }
	  }
	}
      }      
    }

    {
      Gnuplot plot_obj;

      FUNC_LIB::function<scalar_type, w_dmn_t> phi_w("phi_w");

      for(int m2=0; m2<b::dmn_size(); m2++){
	for(int m1=0; m1<b::dmn_size(); m1++){
	  for(int n2=0; n2<b::dmn_size(); n2++){
	    for(int n1=0; n1<b::dmn_size(); n1++){
	      for(int k_ind=0; k_ind<K_dmn::dmn_size(); k_ind++){
		
		if(abs(K_dmn::get_elements()[k_ind][0]+K_dmn::get_elements()[k_ind][1]-M_PI)<1.e-6){
		  for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++)
		    phi_w(w_ind) = imag(chi(n1,n2,m1,m2,k_ind,w_ind));

		  SHOW::execute(plot_obj, phi_w);
		}
	      }
	    }
	  }
	}
      }      
    }

    {
      std::vector<double> x(0);
      std::vector<double> y(0);
      std::vector<double> z(0);

      for(int k_ind=0; k_ind<K_dmn::dmn_size(); k_ind++){
	x.push_back(K_dmn::get_elements()[k_ind][0]);
	y.push_back(K_dmn::get_elements()[k_ind][1]);
	z.push_back(real(phi(0,0,0,0,k_ind,w_dmn_t::dmn_size()/2)));
      }

      SHOW::heatmap(x,y,z);
    }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename w_dmn_t>
  void coarsegraining_tp<parameters_type, K_dmn>::compute_tp(FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu , nu , k_HOST         > >& H_k,
                                                               FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu , nu , k_HOST, w      > >& S_k_w,
                                                               FUNC_LIB::function<std::complex<scalar_type>, dmn_4<b_b, b_b, K_dmn, w_dmn_t> >& chi)
  {
    chi = 0.;

    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_HOST   > > A_k  ("A_k");
    FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, k_HOST, w> > A_k_w("A_k_w");

    transform_to_alpha::forward(1., S_k_w, A_k_w);

    K_dmn               k_domain;
    std::pair<int, int> bounds = concurrency.get_bounds(k_domain);

    for(int k_ind=bounds.first; k_ind<bounds.second; k_ind++)
      {
        for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++)
          {
            int w_1, w_2;

            find_w1_and_w2(w_dmn_t::get_elements(), w_ind, w_1, w_2);

            {
              for(int k=0; k<k_HOST::dmn_size(); k++)
                for(int j=0; j<nu::dmn_size(); j++)
                  for(int i=0; i<nu::dmn_size(); i++)
                    A_k(i,j,k) = A_k_w(i,j,k,w_1);

              coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(k_ind, w_1, H_k, A_k, I_q, H_q, A_q, S_q, G_q);
            }

            {
              for(int k=0; k<k_HOST::dmn_size(); k++)
                for(int j=0; j<nu::dmn_size(); j++)
                  for(int i=0; i<nu::dmn_size(); i++)
                    A_k(i,j,k) = A_k_w(i,j,k,w_2);

              coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(k_ind, w_2, H_k, A_k, I_q_plus_Q, H_q_plus_Q, A_q_plus_Q, S_q_plus_Q, G_q_plus_Q);
            }

            compute_bubble(bubble_q);

            {
              double factor = get_integration_factor();

              for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
                for(int n1=0; n1<b::dmn_size(); n1++)
                  for(int n2=0; n2<b::dmn_size(); n2++)
                    for(int m1=0; m1<b::dmn_size(); m1++)
                      for(int m2=0; m2<b::dmn_size(); m2++)
                        chi(n1,n2,m1,m2,k_ind,w_ind) += factor*w_q(q_ind)*bubble_q(n1,n2,m1,m2,q_ind);
            }
          }
      }

    concurrency.sum(chi);

    {
      scalar_type V_K=0;
      for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
        V_K += w_q(q_ind);

      chi /= V_K;
    }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename w_dmn_t>
  void coarsegraining_tp<parameters_type, K_dmn>::compute_phi(FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu , nu , k_HOST         > >& H_k,
                                                               FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu , nu , k_HOST, w      > >& S_k_w,
                                                               FUNC_LIB::function<std::complex<scalar_type>, dmn_4<b_b, b_b, K_dmn , w_dmn_t> >& phi)
  {
    if(concurrency.id()==concurrency.last())
      cout << "\n\n\t start " << __FUNCTION__ << " ... " << print_time();

    phi = 0.;

    FUNC_LIB::function<std::complex<scalar_type>, dmn_3<nu, nu, k_HOST   > > A_k  ("A_k");
    FUNC_LIB::function<std::complex<scalar_type>, dmn_4<nu, nu, k_HOST, w> > A_k_w("A_k_w");

    transform_to_alpha::forward(1., S_k_w, A_k_w);

    K_dmn               k_domain;
    std::pair<int, int> bounds = concurrency.get_bounds(k_domain);

    for(int k_ind=bounds.first; k_ind<bounds.second; k_ind++)
      {
        for(int w_ind=0; w_ind<w_dmn_t::dmn_size(); w_ind++)
          {
            int w_1, w_2;

            find_w1_and_w2(w_dmn_t::get_elements(), w_ind, w_1, w_2);

            {
              for(int k=0; k<k_HOST::dmn_size(); k++)
                for(int j=0; j<nu::dmn_size(); j++)
                  for(int i=0; i<nu::dmn_size(); i++)
                    A_k(i,j,k) = A_k_w(i,j,k,w_1);

              coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(k_ind, w_1, H_k, A_k, I_q, H_q, A_q, S_q, G_q);
            }

            {
              for(int k=0; k<k_HOST::dmn_size(); k++)
                for(int j=0; j<nu::dmn_size(); j++)
                  for(int i=0; i<nu::dmn_size(); i++)
                    A_k(i,j,k) = A_k_w(i,j,k,w_2);

              coarsegraining_routines<parameters_type, K_dmn>::compute_G_q_w(k_ind, w_2, H_k, A_k, I_Q_min_q, H_Q_min_q, A_Q_min_q, S_Q_min_q, G_Q_min_q);
            }

            compute_bubble(bubble_q);

            {
              double factor = get_integration_factor();

              for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
		for(int m2=0; m2<b::dmn_size(); m2++)
		  for(int m1=0; m1<b::dmn_size(); m1++)
		    for(int n2=0; n2<b::dmn_size(); n2++)
		      for(int n1=0; n1<b::dmn_size(); n1++)
                        phi(n1,n2,m1,m2,k_ind,w_ind) += factor*w_q(q_ind)*bubble_q(n1,n2,m1,m2,q_ind);
            }
          }
      }

    concurrency.sum(phi);

    {
      scalar_type V_K=0;
      for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
        V_K += w_q(q_ind);

      phi /= V_K;
    }

    if(concurrency.id()==concurrency.last())
      cout << "\n\n\t end  ... " << print_time();
  }

  template<typename parameters_type, typename K_dmn>
  void coarsegraining_tp<parameters_type, K_dmn>::find_w1_and_w2(std::vector<double>& elements, int& w_ind, int& w1, int& w2)
  {
    //int N_w   = elements.size();
    int W_ind = parameters.get_w_channel();

    for(int l=0; l<w::dmn_size(); l++)
      if(abs(elements[w_ind]-w::get_elements()[l])<1.e-6)
        w1 = l;

    assert(abs(w::get_elements()[w1]-elements[w_ind]) < 1.e-6);

    switch(parameters.get_vertex_measurement_type())
      {
      case PARTICLE_HOLE_CHARGE:
      case PARTICLE_HOLE_MAGNETIC:
      case PARTICLE_HOLE_TRANSVERSE:
        {
          w2 = w1 + W_ind;
          //assert(abs(w::get_elements()[w2] - (elements[N_w/2+W]+elements[w_ind])) < 1.e-6);
        }
        break;

      case PARTICLE_PARTICLE_SUPERCONDUCTING:
        {
          w2 = W_ind + (w::dmn_size()-1-w1);
          assert(abs(w::get_elements()[w1]+w::get_elements()[w::dmn_size()-1-w1]) < 1.e-6);
          //assert(abs(w::get_elements()[w2] - (elements[N_w/2+W]-elements[w_ind]+)) < 1.e-6);
        }
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }

  template<typename parameters_type, typename K_dmn>
  void coarsegraining_tp<parameters_type, K_dmn>::compute_bubble(FUNC_LIB::function<std::complex<scalar_type>, dmn_3<b_b, b_b, q_dmn> >& bubble)
  {
    bubble = 0.;

    for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++){
      for(int n1=0; n1<b::dmn_size(); n1++){
        for(int n2=0; n2<b::dmn_size(); n2++){
          for(int m1=0; m1<b::dmn_size(); m1++){
            for(int m2=0; m2<b::dmn_size(); m2++){

              switch(parameters.get_vertex_measurement_type())
                {
                case PARTICLE_HOLE_TRANSVERSE:
                  bubble(n1,n2,m1,m2,q_ind) += G_q(n1, e_UP, m2, e_UP, q_ind) * G_q_plus_Q(n2, e_UP, m1, e_UP, q_ind);
                  break;

                case PARTICLE_HOLE_MAGNETIC:
                  bubble(n1,n2,m1,m2,q_ind) += G_q(n1, e_UP, m2, e_UP, q_ind) * G_q_plus_Q(n2, e_UP, m1, e_UP, q_ind);
                  break;

                case PARTICLE_HOLE_CHARGE:
                  bubble(n1,n2,m1,m2,q_ind) += G_q(n1, e_UP, m2, e_UP, q_ind) * G_q_plus_Q(n2, e_UP, m1, e_UP, q_ind);
                  break;

                case PARTICLE_PARTICLE_SUPERCONDUCTING:
                  bubble(n1,n2,m1,m2,q_ind) += G_q(n1, e_UP, m1, e_UP, q_ind) * G_Q_min_q(n2, e_UP, m2, e_UP, q_ind);
                  break;

                default:
                  throw std::logic_error(__FUNCTION__);
                }
            }
          }
        }
      }
    }
  }

  template<typename parameters_type, typename K_dmn>
  double coarsegraining_tp<parameters_type, K_dmn>::get_integration_factor()
  {
    switch(parameters.get_vertex_measurement_type())
      {
      case PARTICLE_HOLE_TRANSVERSE:
        return -1.;
        break;

      case PARTICLE_HOLE_MAGNETIC:
        return -1.;
        break;

      case PARTICLE_HOLE_CHARGE:
        return -2.;
        break;

      case PARTICLE_PARTICLE_SUPERCONDUCTING:
        return 1;
        break;

      default:
        throw std::logic_error(__FUNCTION__);
      }
  }

  /*
    template<typename parameters_type, typename K_dmn>
    std::pair<int, std::vector<double> > coarsegraining_tp<parameters_type, K_dmn>::get_k_accent(int K_ind)
    {
    //int Q_ind             = parameters.get_q_channel();
    std::vector<double> Q = parameters.get_q_vector();

    int Q_ind=-1;
    for(int l=0; l<k_cluster_type::get_size(); ++l)
    if(VECTOR_OPERATIONS::L2_NORM(Q, k_cluster_type::get_elements()[l])<1.e-6)
    Q_ind = l;

    assert(Q_ind>-1);

    std::pair<int, std::vector<double> > k_accent;

    switch(parameters.get_vertex_measurement_type())
    {
    case PARTICLE_HOLE_TRANSVERSE:
    {
    int index                    = k_cluster_type::add(K_ind, Q_ind);
    std::vector<double> K_plus_Q = k_cluster_type::get_elements()[K_ind];

    for(int i=0; i<int(K_plus_Q.size()); i++)
    K_plus_Q[i] += Q[i];

    k_accent.first  = index;
    k_accent.second = K_plus_Q;
    }
    break;

    case PARTICLE_HOLE_MAGNETIC:
    {
    int index                    = k_cluster_type::add(K_ind, Q_ind);
    std::vector<double> K_plus_Q = k_cluster_type::get_elements()[K_ind];

    for(int i=0; i<int(K_plus_Q.size()); i++)
    K_plus_Q[i] += Q[i];

    k_accent.first  = index;
    k_accent.second = K_plus_Q;
    }
    break;

    case PARTICLE_HOLE_CHARGE:
    {
    int index                    = k_cluster_type::add(K_ind, Q_ind);
    std::vector<double> K_plus_Q = k_cluster_type::get_elements()[K_ind];

    for(int i=0; i<int(K_plus_Q.size()); i++)
    K_plus_Q[i] += Q[i];

    k_accent.first  = index;
    k_accent.second = K_plus_Q;
    }
    break;

    case PARTICLE_PARTICLE_SUPERCONDUCTING:
    {
    int index                   = k_cluster_type::subtract(K_ind, Q_ind);
    std::vector<double> Q_min_K = Q;

    for(int i=0; i<int(Q_min_K.size()); i++)
    Q_min_K[i] -= k_cluster_type::get_elements()[K_ind][i];

    k_accent.first  = index;
    k_accent.second = Q_min_K;
    }
    break;

    default:
    throw std::logic_error(__FUNCTION__);
    }

    return k_accent;
    }

    template<typename parameters_type, typename K_dmn>
    std::pair<int, double> coarsegraining_tp<parameters_type, K_dmn>::get_w_accent(int w_index, int W)
    {
    std::pair<int, double> w_accent;

    switch(parameters.get_vertex_measurement_type())
    {
    case PARTICLE_HOLE_TRANSVERSE:
    w_accent.first = w_index + W;
    assert(w_accent.first >= 0 && w_accent.first < frequency_domain_type::get_size());
    w_accent.second = frequency_domain_type::get_elements()[w_accent.first];
    break;

    case PARTICLE_HOLE_MAGNETIC:
    w_accent.first = w_index + W;
    assert(w_accent.first >= 0 && w_accent.first < frequency_domain_type::get_size());
    w_accent.second = frequency_domain_type::get_elements()[w_accent.first];
    break;

    case PARTICLE_HOLE_CHARGE:
    w_accent.first = w_index + W;
    assert(w_accent.first >= 0 && w_accent.first < frequency_domain_type::get_size());
    w_accent.second = frequency_domain_type::get_elements()[w_accent.first];
    break;

    case PARTICLE_PARTICLE_SUPERCONDUCTING:
    w_accent.first = W + (frequency_domain_type::get_size()-1)-w_index ;
    assert(w_accent.first >= 0 && w_accent.first < frequency_domain_type::get_size());
    w_accent.second = frequency_domain_type::get_elements()[w_accent.first];
    break;

    default:
    throw std::logic_error(__FUNCTION__);
    }

    return w_accent;
    }
  */

}



#endif
