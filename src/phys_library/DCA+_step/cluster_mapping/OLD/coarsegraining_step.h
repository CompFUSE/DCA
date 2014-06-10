//-*-C++-*-

#ifndef DCA_CLUSTER_MAP_H
#define DCA_CLUSTER_MAP_H

namespace DCA
{
//   enum COARSEGRAIN_DOMAIN_NAMES {ORIGIN, K, K_PLUS_Q, Q_MINUS_K};

//   std::string to_str(COARSEGRAIN_DOMAIN_NAMES NAME)
//   {
//     switch(NAME)
//       {
//       case ORIGIN:
// 	return "ORIGIN";

//       case K:
// 	return "K";

//       case K_PLUS_Q:
// 	return "K_PLUS_Q";

//       case Q_MINUS_K:
// 	return "Q_minus_K";

//       default:
// 	throw std::logic_error(__FUNCTION__);
//       }

//     return "???";
//   }
  
//   template<typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
//   class coarsegrain_domain
//   {  
//   public:
    
//     const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

//     typedef double              scalar_type;
//     typedef std::vector<double> element_type;
    
//     typedef MATH_ALGORITHMS::domain_specifications<scalar_type, element_type, 
// 						   MATH_ALGORITHMS::DISCRETE, MATH_ALGORITHMS::KRONECKER_DELTA, 
// 						   MATH_ALGORITHMS::INTERVAL, MATH_ALGORITHMS::EQUIDISTANT>     dmn_specifications_type;

//   public:
    
//     static int& get_size()
//     {
//       static int size = 0;
//       return size;
//     }
    
//     static std::string get_name()
//     {
//       static std::string name = "coarsegrain_domain (" + to_str(NAME) + ")";
//       return name;
//     }

//     static std::vector<scalar_type>& get_weights()
//     {
//       static std::vector<scalar_type> weights(0);
//       return weights;
//     }

//     static std::vector<element_type>& get_elements()
//     {
//       static std::vector<element_type> elements(0);
//       return elements;
//     }

//     static void set_elements(int K_ind)
//     {
//       std::vector<element_type> elements = coarsegrain_domain<K_dmn, ORIGIN>::get_elements();

//       switch(NAME)
// 	{
// 	case K:
// 	  {
// 	    for(int q_ind=0; q_ind<elements.size(); q_ind++)
// 	      for(int d_ind=0; d_ind<DIMENSION; d_ind++)
// 		elements[q_ind][d_ind] += K_dmn::get_elements()[K_ind][d_ind];
// 	  }
// 	  break;

// 	default:
// 	  throw std::logic_error(__FUNCTION__);
// 	}

//       get_elements() = elements;
//     }

//     static std::vector<element_type>& set_elements(int K_ind, int Q_ind)
//     {
//       std::vector<element_type> elements = coarsegrain_domain<K_dmn, ORIGIN>::get_elements();

//       switch(NAME)
// 	{
// 	case K_PLUS_Q:
// 	  {
// 	    for(int q_ind=0; q_ind<elements.size(); q_ind++){
// 	      for(int d_ind=0; d_ind<DIMENSION; d_ind++){
// 		elements[q_ind][d_ind] += K_dmn::get_elements()[K_ind][d_ind];
// 		elements[q_ind][d_ind] += K_dmn::get_elements()[Q_ind][d_ind];
// 	      }
// 	    }
// 	  }
// 	  break;

// 	case Q_MINUS_K:
// 	  {
// 	    for(int q_ind=0; q_ind<elements.size(); q_ind++){
// 	      for(int d_ind=0; d_ind<DIMENSION; d_ind++){
// 		elements[q_ind][d_ind] *= -1;
// 		elements[q_ind][d_ind] -= K_dmn::get_elements()[K_ind][d_ind];
// 		elements[q_ind][d_ind] += K_dmn::get_elements()[Q_ind][d_ind];
// 	      }
// 	    }
// 	  }
// 	  break;

// 	default:
// 	  throw std::logic_error(__FUNCTION__);
// 	}

//       get_elements() = elements;
//     }
    
//   };

  template<typename k_dmn, typename K_dmn>
  class interpolation_matrices
  {};

  template<typename k_dmn, typename K_dmn, COARSEGRAIN_DOMAIN_NAMES NAME>
  class interpolation_matrices<k_dmn, dmn_0<coarsegrain_domain<K_dmn, NAME> > >
  {  
    typedef LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> matrix_type;

    typedef dmn_0<coarsegrain_domain<K_dmn, NAME> >             q_dmn;
    typedef typename k_dmn::parameter_type::dual_type           r_dmn;

    typedef MATH_ALGORITHMS::basis_transform<typename k_dmn::parameter_type, r_dmn> trafo_k_to_r_type;
    typedef MATH_ALGORITHMS::basis_transform<r_dmn, typename q_dmn::parameter_type> trafo_r_to_q_type;

  public:
    
    static function<matrix_type, K_dmn>& get();

    static matrix_type& get(int k_ind);

    static bool is_initialized();

    template<typename parameters_type>
    static void initialize(parameters_type& parameters);
  };

  template<typename k_dmn, typename K_dmn>
  class interpolation_matrices<k_dmn, dmn_0<coarsegrain_domain<K_dmn, K> > >
  {  
    typedef LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU> matrix_type;

    typedef dmn_0<coarsegrain_domain<K_dmn, K> >      q_dmn;
    typedef typename k_dmn::parameter_type::dual_type r_dmn;

    typedef MATH_ALGORITHMS::basis_transform<typename k_dmn::parameter_type, r_dmn> trafo_k_to_r_type;
    typedef MATH_ALGORITHMS::basis_transform<r_dmn, typename q_dmn::parameter_type> trafo_r_to_q_type;

  public:
    
    static function<matrix_type, K_dmn>& get()
    {
      static function<matrix_type, K_dmn> k_to_q("k_to_q ("+q_dmn::parameter_type::get_name()+")");
      assert(is_initialized()==true);

      return k_to_q;
    }

    static matrix_type& get(int k_ind)
    {      
      static function<matrix_type, K_dmn>& k_to_q = get();
      assert(is_initialized()==true);

      return k_to_q(k_ind);
    }

    static bool is_initialized()
    {
      static bool initialized = false;
      return initialized;
    }

    template<typename parameters_type>
    static void initialize(parameters_type& parameters);
  };

  template<typename k_dmn, typename K_dmn>
  template<typename parameters_type>
  void interpolation_matrices<k_dmn, dmn_0<coarsegrain_domain<K_dmn, K> > >::initialize(parameters_type& parameters)
  {
    is_initialized() = true;

    typedef typename parameters_type::Concurrency_Type concurrency_type;

    concurrency_type& concurrency = parameters.get_concurrency();

    K_dmn K_dmn_obj;
    std::pair<int, int> bounds = concurrency.get_bounds(K_dmn_obj);
      
    for(int K_ind=bounds.first; K_ind<bounds.second; K_ind++)
      {
	q_dmn::parameter_type::set_elements(K_ind);

	matrix_type& T_k_to_q = get(K_ind);

	T_k_to_q.resize(std::pair<int, int>(q_dmn::dmn_size(), k_dmn::dmn_size()));

	trafo_k_to_r_type::is_initialized() = false;
	trafo_r_to_q_type::is_initialized() = false;

	matrix_type& T_k_to_r = trafo_k_to_r_type::get_transformation_matrix();  
	matrix_type& T_r_to_q = trafo_r_to_q_type::get_transformation_matrix();  

	LIN_ALG::GEMM<LIN_ALG::CPU>::execute(T_r_to_q, T_k_to_r, T_k_to_q);
      }
    
    for(int K_ind=0; K_ind<K_dmn::dmn_size(); K_ind++)
      concurrency.sum(get(K_ind));
  }


  template<typename parameters_type, typename K_dmn>
  class cluster_map
  {
#include "type_definitions.h"
    
    const static int DIMENSION = K_dmn::parameter_type::DIMENSION;

    typedef typename parameters_type::Concurrency_Type concurrency_type;

    typedef double scalar_type;

    typedef dmn_0<tetrahedron_mesh<K_dmn> >             tetrahedron_dmn;

    typedef MATH_ALGORITHMS::gaussian_quadrature_domain<tetrahedron_dmn> quadrature_dmn;

    typedef dmn_0<coarsegrain_domain<K_dmn, ORIGIN   > > q_0_dmn;

    typedef dmn_0<coarsegrain_domain<K_dmn, K        > > q_dmn;
    typedef dmn_0<coarsegrain_domain<K_dmn, K_PLUS_Q > > q_plus_k_dmn;
    typedef dmn_0<coarsegrain_domain<K_dmn, Q_MINUS_K> > k_minus_q_dmn;

    typedef dmn_3<nu, nu, q_dmn>                         nu_nu_q;

    typedef typename parameters_type::parallelization_t parallelization_type;

    typedef LIN_ALG::matrix<double, LIN_ALG::CPU>       matrix_type;

  public:

    cluster_map(parameters_type& parameters_ref);
  
    ~cluster_map();

    void initialize();

    template<typename k_dmn>
    void plot_H_q(function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn> >& H_0);

    template<typename k_dmn>
    void plot_S_q(function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w);
		
    template<typename k_dmn>
    void compute_S_K_w(function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w,
		       function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn, w> >& S_K_w);

    void compute_G0_K_t(function<std::complex<scalar_type>, dmn_3<nu, nu, k_LDA   > >& H_0,
			function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn, t> >& G0_k_w);

    template<typename k_dmn>
    void compute_G_K_w(function<std::complex<scalar_type>, dmn_3<nu, nu, k_LDA> >&    H_0,
		       function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w,
		       function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn, w> >& G_K_w);
		      
  private:
  
    void compute_gaussian_mesh();

    template<typename tmp_scalar_type>
    void compute_I_q(tmp_scalar_type value);
    
    template<typename k_dmn>
    void compute_H_q  (int K_ind,            function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn> >& H_0);
    
    template<typename k_dmn>
    void compute_S_q_w(int K_ind, int w_ind, function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w);
    
    void compute_S_q_w(int K_ind, int w_ind, function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn, w> >& S_k_w);
    
    template<typename k_dmn>
    void compute_G_q_t(int K_ind, int t_ind, 
		       function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn> >& H_0);

    template<typename k_dmn>
    void compute_G_q_w(int K_ind, int w_ind, 
		       function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn   > >& H_0,
		       function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w);
    
  private:

    parameters_type&  parameters;
    concurrency_type& concurrency;
  
    //   matrix_type k_to_K;

    //   function<matrix_type, K_dmn> k_to_q;

    //   function<std::complex<scalar_type>, nu_nu_k> S_K;
    //   function<std::complex<scalar_type>, nu_nu_k> G_K;

    //   function<std::complex<scalar_type>, nu_nu_k> H_k;
    //   function<std::complex<scalar_type>, nu_nu_k> S_k;
    //   function<std::complex<scalar_type>, nu_nu_k> G_k;

    function<std::complex<scalar_type>, nu_nu_q> I_q;
    function<std::complex<scalar_type>, nu_nu_q> H_q;
    function<std::complex<scalar_type>, nu_nu_q> S_q;
    function<std::complex<scalar_type>, nu_nu_q> G_q;
  };

  template<typename parameters_type, typename K_dmn>
  cluster_map<parameters_type, K_dmn>::cluster_map(parameters_type& parameters_ref):
    parameters (parameters_ref),
    concurrency(parameters.get_concurrency())// ,

    //   k_to_K("k_to_K", std::pair<int, int>(K_dmn::dmn_size(), k_dmn::dmn_size()))
  {}

  template<typename parameters_type, typename K_dmn>
  cluster_map<parameters_type, K_dmn>::~cluster_map()
  {}

  template<typename parameters_type, typename K_dmn>
  void cluster_map<parameters_type, K_dmn>::initialize()
  {
    //SHOW::plot_points(K_dmn::get_elements());

    compute_gaussian_mesh();

    //compute_interpolation_matrices_with_gaussian_fitting();

    {
      I_q.reset();
      H_q.reset();
      S_q.reset();
      G_q.reset();
    }
  }

  // template<typename parameters_type, typename K_dmn>
  // double cluster_map<parameters_type, K_dmn>::compute_length_scale()
  // {
  //   double L=0;

  //   {
  //     std::vector<double> L_d(DIMENSION, 0); 

  //     for(int j=0; j<DIMENSION; j++)
  //       L_d[j] = 0;
      
  //     for(int j=0; j<DIMENSION; j++)
  //       for(int i=0; i<DIMENSION; i++)
  // 	L_d[j] += std::pow(k_dmn::parameter_type::get_basis()[i+j*DIMENSION], 2);
      
  //     for(int j=0; j<DIMENSION; j++)
  //       L_d[j] = sqrt(L_d[j]);

  //     for(int i=0; i<DIMENSION; i++)
  //       L += L_d[i]/scalar_type(DIMENSION);

  //     VECTOR_OPERATIONS::PRINT(L_d);
  //     cout << endl;
  //   }

  //   return L;
  // }

  template<typename parameters_type, typename K_dmn>
  void cluster_map<parameters_type, K_dmn>::compute_gaussian_mesh()
  {
    quadrature_dmn::initialize_Brillouin_zone(parameters.get_mesh_refinement(),
					      parameters.get_quadrature_rule(),
					      parameters.get_phi_period());
    
    //SHOW::plot_points(quadrature_dmn::get_elements());
  
    {
      q_0_dmn::parameter_type::get_size()     = quadrature_dmn::get_size();
      q_0_dmn::parameter_type::get_weights()  = quadrature_dmn::get_weights();
      q_0_dmn::parameter_type::get_elements() = quadrature_dmn::get_elements();
    }

    {
      q_dmn::parameter_type::get_size()     = quadrature_dmn::get_size();    
      q_dmn::parameter_type::get_weights()  = quadrature_dmn::get_weights();    
      q_dmn::parameter_type::get_elements() = quadrature_dmn::get_elements();
    }

    {
      q_plus_k_dmn::parameter_type::get_size()     = quadrature_dmn::get_size();    
      q_plus_k_dmn::parameter_type::get_weights()  = quadrature_dmn::get_weights();    
      q_plus_k_dmn::parameter_type::get_elements() = quadrature_dmn::get_elements();
    }

    {
      k_minus_q_dmn::parameter_type::get_size()     = quadrature_dmn::get_size();    
      k_minus_q_dmn::parameter_type::get_weights()  = quadrature_dmn::get_weights();    
      k_minus_q_dmn::parameter_type::get_elements() = quadrature_dmn::get_elements();
    }
  }

  /*
    template<typename parameters_type, typename K_dmn>
    void cluster_map<parameters_type, K_dmn>::plot_H_q(function<std::complex<scalar_type>, nu_nu_k>& H_0)
    {
    function<double, k_dmn> tmp_0;
    function<double, q_dmn> tmp_1;

    std::vector<double> x(q_dmn::dmn_size()*K_dmn::dmn_size());
    std::vector<double> y(q_dmn::dmn_size()*K_dmn::dmn_size());
    std::vector<double> z(q_dmn::dmn_size()*K_dmn::dmn_size());

    gaussian_fit<scalar_type, k_dmn, q_dmn> gaussian_fit_k_to_q;

    for(int K_ind=0; K_ind<K_dmn::dmn_size(); K_ind++){

    LIN_ALG::matrix<double, LIN_ALG::CPU>& T = k_to_q(K_ind);
    
    double* A_ptr = &real(H_0(0));
    double* B_ptr = &T(0,0);
    double* C_ptr = &real(H_q(0));
    
    int M = 2*nu_nu::dmn_size();
    int K = k_dmn::dmn_size();
    int N = q_dmn::dmn_size();
    
    int LDA = 2*nu_nu::dmn_size();
    int LDB = T.get_global_size().first;
    int LDC = 2*nu_nu::dmn_size();
    
    LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'T', M, N, K, 1., A_ptr, LDA, B_ptr, LDB, 0., C_ptr, LDC);
    
    for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++){
    x[q_ind+K_ind*q_dmn::dmn_size()] = q_0_dmn::parameter_type::get_elements()[q_ind][0]+K_dmn::get_elements()[K_ind][0];
    y[q_ind+K_ind*q_dmn::dmn_size()] = q_0_dmn::parameter_type::get_elements()[q_ind][1]+K_dmn::get_elements()[K_ind][1];
    z[q_ind+K_ind*q_dmn::dmn_size()] = real(H_q(0,0,q_ind));
    }
    }
  
    //   SHOW::plot_points(x,y);
    //   SHOW::heatmap(x,y,z);
    }
  */

  template<typename parameters_type, typename K_dmn>
  template<typename tmp_scalar_type>
  void cluster_map<parameters_type, K_dmn>::compute_I_q(tmp_scalar_type value)
  { 
    for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
      for(int j=0; j<nu::dmn_size(); j++)
	for(int i=0; i<nu::dmn_size(); i++)
	  I_q(i,j,q_ind) = i==j? value : 0.;
  }

  template<typename parameters_type, typename K_dmn>
  template<typename k_dmn>
  void cluster_map<parameters_type, K_dmn>::compute_H_q(int K_ind, function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn> >& H_0)
  {
    //   LIN_ALG::matrix<double, LIN_ALG::CPU>& T = k_to_q(K_ind);
    
    //   double* A_ptr = &real(H_0(0));
    //   double* B_ptr = &T(0,0);
    //   double* C_ptr = &real(H_q(0));
  
    //   int M = 2*nu_nu::dmn_size();
    //   int K = k_dmn::dmn_size();
    //   int N = q_dmn::dmn_size();
  
    //   int LDA = 2*nu_nu::dmn_size();
    //   int LDB = T.get_global_size().first;
    //   int LDC = 2*nu_nu::dmn_size();
  
    //   LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'T', M, N, K, 1., A_ptr, LDA, B_ptr, LDB, 0., C_ptr, LDC);
  }

  template<typename parameters_type, typename K_dmn>
  template<typename k_dmn>
  void cluster_map<parameters_type, K_dmn>::compute_S_q_w(int K_ind, int w_ind, function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w)
  {
    //   LIN_ALG::matrix<double, LIN_ALG::CPU>& T = k_to_q(K_ind);
    
    //   double* A_ptr = &real(S_k_w(0, 0, 0, w_ind));
    //   double* B_ptr = &T(0,0);
    //   double* C_ptr = &real(S_q(0));
  
    //   int M = 2*nu_nu::dmn_size();
    //   int K = k_dmn::dmn_size();
    //   int N = q_dmn::dmn_size();
  
    //   int LDA = 2*nu_nu::dmn_size();
    //   int LDB = T.get_global_size().first;
    //   int LDC = 2*nu_nu::dmn_size();
  
    //   LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'T', M, N, K, 1., A_ptr, LDA, B_ptr, LDB, 0., C_ptr, LDC);
  }

  template<typename parameters_type, typename K_dmn>
  void cluster_map<parameters_type, K_dmn>::compute_S_q_w(int K_ind, int w_ind, function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn, w> >& S_K_w)
  {
    for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
      memcpy(&S_q(0,0,q_ind), &S_K_w(0, 0, K_ind, w_ind), sizeof(std::complex<scalar_type>)*nu_nu::dmn_size());
  }

  template<typename parameters_type, typename K_dmn>
  template<typename k_dmn>
  void cluster_map<parameters_type, K_dmn>::compute_G_q_w(int K_ind, int w_ind, 
							  function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn   > >&    H_0,
							  function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w)
  {
    {
      std::complex<double> I(0, 1); 
      std::complex<double> i_wm_min_mu = w::get_elements()[w_ind]*I+parameters.get_chemical_potential();

      compute_I_q(i_wm_min_mu);
    }

    compute_H_q(K_ind, H_0);

    compute_S_q_w(K_ind, w_ind, S_k_w);

    LIN_ALG::matrix<std::complex<scalar_type>, LIN_ALG::CPU> G_inv("G_inv", nu::dmn_size());

    for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++){
    
      for(int j=0; j<nu::dmn_size(); j++)
	for(int i=0; i<nu::dmn_size(); i++)
	  G_inv(i,j) = I_q(i,j,q_ind)-H_q(i,j,q_ind)-S_q(i,j,q_ind);

      LIN_ALG::GEINV<LIN_ALG::CPU>::execute(G_inv);

      for(int j=0; j<nu::dmn_size(); j++)
	for(int i=0; i<nu::dmn_size(); i++)
	  G_q(i,j,q_ind) = G_inv(i,j);
    }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename k_dmn>
  void cluster_map<parameters_type, K_dmn>::compute_G_K_w(function<std::complex<scalar_type>, dmn_3<nu, nu, k_LDA> >&    H_0,
							  function<std::complex<scalar_type>, dmn_4<nu, nu, k_dmn, w> >& S_k_w,
							  function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn, w> >& G_K_w)
  {
    cout << __FUNCTION__ << endl;

    G_K_w = 0.;

    std::vector<scalar_type>& weights = q_0_dmn::parameter_type::get_weights();

    dmn_2<K_dmn, w> K_wm_dmn;  
    std::pair<int, int> bounds = concurrency.get_bounds(K_wm_dmn);
  
    int coor[2];
    for(int l=bounds.first; l<bounds.second; l++)
      {
	K_wm_dmn.linind_2_subind(l, coor);

	compute_G_q_w(coor[0], coor[1], H_0, S_k_w);

	for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
	  for(int j=0; j<nu::dmn_size(); j++)
	    for(int i=0; i<nu::dmn_size(); i++)
	      G_K_w(i,j,coor[0],coor[1]) += G_q(i,j,q_ind)*weights[q_ind];
      }

    concurrency.sum(G_K_w);

    {
      double V_K=0;
      for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
	V_K += weights[q_ind];

      G_K_w /= V_K;
    }
  }

  template<typename parameters_type, typename K_dmn>
  template<typename k_dmn>
  void cluster_map<parameters_type, K_dmn>::compute_G_q_t(int K_ind, int t_ind, 
							  function<std::complex<scalar_type>, dmn_3<nu, nu, k_dmn> >& H_0)
  {
    double t_val = t::get_elements()[t_ind]; 
    double beta  = parameters.get_beta();

    compute_I_q(parameters.get_chemical_potential());

    compute_H_q(K_ind, H_0);

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
	G_t[i] = std::exp(L[i]*(beta-t_val))/(std::exp(L[i]*beta)+1.);

      for(int j=0; j<nu::dmn_size(); j++)
	for(int i=0; i<nu::dmn_size(); i++)
	  for(int l=0; l<nu::dmn_size(); l++)
	    G_q(i,j,q_ind) += G_t[l]*real(conj(V(l,i))*V(l,j));
    }
  }

  template<typename parameters_type, typename K_dmn>
  void cluster_map<parameters_type, K_dmn>::compute_G0_K_t(function<std::complex<scalar_type>, dmn_3<nu, nu, k_LDA   > >&   H_0,
							   function<std::complex<scalar_type>, dmn_4<nu, nu, K_dmn, t> >& G_K_t)
							   {
							     cout << __FUNCTION__ << endl;

							     G_K_t = 0.;

							     std::vector<scalar_type>& weights = q_0_dmn::parameter_type::get_weights();

							     dmn_2<K_dmn, t> K_t_dmn;  
							     std::pair<int, int> bounds = concurrency.get_bounds(K_t_dmn);
  
							     int coor[2];
							     for(int l=bounds.first; l<bounds.second; l++)
							       {
								 K_t_dmn.linind_2_subind(l, coor);

								 compute_G_q_t(coor[0], coor[1], H_0);
      
								 for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
								   for(int j=0; j<nu::dmn_size(); j++)
								     for(int i=0; i<nu::dmn_size(); i++)
								       G_K_t(i,j,coor[0],coor[1]) += G_q(i,j,q_ind)*weights[q_ind];
							       }

							     concurrency.sum(G_K_t);

							     {
							       double V_K=0;
							       for(int q_ind=0; q_ind<q_dmn::dmn_size(); q_ind++)
								 V_K += weights[q_ind];

							       G_K_t /= V_K;
							     }
							   }

							   }
							 
#endif
