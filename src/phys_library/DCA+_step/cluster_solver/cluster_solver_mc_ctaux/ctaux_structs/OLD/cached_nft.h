//-*-C++-*-

/** \ingroup DCA */

/*@{*/

/*! \file cached_auxilery_field_values.h  
 *
 *  Contains a templated class over the dimension to represent the cached nft-plan
 */

#ifndef CACHED_NFT_H
#define CACHED_NFT_H

namespace QMC {

  template<typename whatever_1, typename whatever_2, typename whatever_3>
  struct triple
  {
    typedef triple<whatever_1, whatever_2, whatever_3> this_type;

    whatever_1 first;
    whatever_2 second;
    whatever_3 third;

    static bool less(const this_type& l, const this_type& r)
    {
      if(l.first == r.first)
	if(l.second == r.second)
	  return (l.third < r.third);
	else
	  return (l.second < r.second);
      else
	return (l.first < r.first);
    }
  };


  template<int dimension>
  class cached_nft
  {
#include "type_definitions.h" 

  public:

    cached_nft(std::vector<double>& matsubara_frequencies_ref);

    ~cached_nft();

    template<typename configuration_t,
	     typename matrix_t,
	     typename function_t>
    void execute(configuration_t& configuration,
		 matrix_t&        M,
		 function_t&      M_r_r_w_w);

  private:

    void initialize();

    template<typename configuration_t>
    void compute_T(configuration_t& configuration);

    template<typename configuration_t>
    int sort_configuration(configuration_t&                  configuration,
			   std::vector<std::pair<int,int> >& p,
			   std::vector<int>&                 start_index,
			   std::vector<int>&                 end_index);

    template<typename configuration_t>
    int sort_configuration(configuration_t&                  configuration,
			   std::vector<triple<int,int,int> >& p,
			   std::vector<int>&                 start_index,
			   std::vector<int>&                 end_index);

  private:

    class trigoniometric_dmn
    {
      const static int INTERPOLATION_NUMBER = 100;

    public:
      static int get_size()
      { 
	return INTERPOLATION_NUMBER; 
      }

      static double get_elements(int n)
      {
	assert(n>=0 && n<get_size());
	return 2*M_PI*double(n)/double(get_size());
      }
    };

    std::vector<double> W;

    FUNC_LIB::function<double, dmn_0<trigoniometric_dmn> > SIN;
    FUNC_LIB::function<double, dmn_0<trigoniometric_dmn> > COS;

    std::complex<double>* T;

    std::complex<double>* M_ij;

    std::complex<double>* T_l;
    std::complex<double>* T_r;

    std::complex<double>* T_l_times_M_ij;
    std::complex<double>* T_l_times_M_ij_times_T_r;
  };

  template<int dimension>
  cached_nft<dimension>::cached_nft(std::vector<double>& matsubara_frequencies_ref):
    W(matsubara_frequencies_ref)
  {
    initialize();
  }

  template<int dimension>
  cached_nft<dimension>::~cached_nft()
  {}

  template<int dimension>
  void cached_nft<dimension>::initialize()
  {
    for(int l=0; l<trigoniometric_dmn::get_size(); l++){
      SIN(l) = sin(trigoniometric_dmn::get_elements(l));
      COS(l) = COS(trigoniometric_dmn::get_elements(l));
    }
  }

  template<int dimension>
  template<typename configuration_t,
	   typename matrix_t,
	   typename function_t>
  void cached_nft<dimension>::execute(configuration_t& configuration,
				      matrix_t&        M,
				      function_t&      M_r_r_w_w)
  {
    //     clock_t begin = clock();

    static int r0_index = DCA_r_cluster_type::get_r_0_index();
    static int N_r      = DCA_cluster_type::get_cluster_size();
    static int N_b      = b::dmn_size();
    static int N_w      = w_VERTEX::dmn_size();
 
    int N_v = configuration.size();

    T = new std::complex<double>[N_w*N_v];

    compute_T(configuration);

    //std::vector<std::pair<int,int> > p(N_v);
    std::vector<triple<int,int,int> > p(N_v);
    
    std::vector<int>                 start_index(N_b*N_r, 0  );
    std::vector<int>                 end_index  (N_b*N_r, N_v);

    int MAX = sort_configuration(configuration, p, start_index, end_index);

    M_ij = new std::complex<double>[MAX*MAX];
    T_l  = new std::complex<double>[N_w*MAX];
    T_r  = new std::complex<double>[N_w*MAX];

    T_l_times_M_ij           = new std::complex<double>[N_w*MAX];
    T_l_times_M_ij_times_T_r = new std::complex<double>[N_w*N_w];

    for(int b_i=0; b_i<N_b; b_i++){
      
      for(int r_i=0; r_i<N_r; r_i++){
	
	int n_I = end_index[b_i + N_b*r_i]-start_index[b_i + N_b*r_i];
	
	for(int b_j=0; b_j<N_b; b_j++){
	  
	  for(int r_j=0; r_j<N_r; r_j++){

	    int min_r_j = DCA_r_cluster_type::subtract(r_j, r0_index);
	    
	    int n_J = end_index[b_j + N_b*r_j]-start_index[b_j + N_b*r_j];
	    
	    if(n_I > 0 && n_J > 0)
	      {

		// M_ij-matrix
		for(int l_i=start_index[b_i + N_b*r_i]; l_i<end_index[b_i + N_b*r_i]; l_i++){

		  assert(p[l_i].first                           == b_i);
		  assert(configuration[p[l_i].third].get_band() == b_i);
	      
		  assert(p[l_i].second                            == r_i);
		  assert(configuration[p[l_i].third].get_r_site() == r_i);

		  int I = l_i-start_index[b_i + N_b*r_i];
	      
		  for(int l_j=start_index[b_j + N_b*r_j]; l_j<end_index[b_j + N_b*r_j]; l_j++){

		    assert(p[l_j].first                           == b_j);
		    assert(configuration[p[l_j].third].get_band() == b_j);
		    
		    assert(p[l_j].second                            == r_j);
		    assert(configuration[p[l_j].third].get_r_site() == r_j);
		    
		    int J = l_j-start_index[b_j + N_b*r_j];
		  
		    M_ij[I + MAX*J] = M( p[l_i].third, p[l_j].third);
		  }
		}

		// T_l matrix
		for(int l_i=start_index[b_i + N_b*r_i]; l_i<end_index[b_i + N_b*r_i]; l_i++){
		  int I = l_i-start_index[b_i + N_b*r_i];
		  memcpy(&T_l[I*N_w], &T[(p[l_i].third)*N_w], sizeof(std::complex<double>)*N_w);
		}

		// T_r matrix
		for(int l_j=start_index[b_j + N_b*r_j]; l_j<end_index[b_j + N_b*r_j]; l_j++){
		  int J = l_j-start_index[b_j + N_b*r_j];
		  memcpy(&T_r[J*N_w], &T[(p[l_j].third)*N_w], sizeof(std::complex<double>)*N_w);
		}
	
		{
		  gemm_plan<std::complex<double> > zgemm_pln;
	  
		  zgemm_pln.A = T_l;
		  zgemm_pln.B = M_ij;
		  zgemm_pln.C = T_l_times_M_ij;

		  zgemm_pln.M = N_w;
		  zgemm_pln.K = n_I;
		  zgemm_pln.N = n_J;

		  zgemm_pln.LDA = N_w;
		  zgemm_pln.LDB = MAX;
		  zgemm_pln.LDC = N_w;

		  zgemm_pln.execute_plan();
		}

		{
		  gemm_plan<std::complex<double> > zgemm_pln;
	  
		  zgemm_pln.TRANSA = 'N';
		  zgemm_pln.TRANSB = 'C';

		  zgemm_pln.A = T_l_times_M_ij;
		  zgemm_pln.B = T_r;
		  zgemm_pln.C = T_l_times_M_ij_times_T_r;

		  zgemm_pln.M = N_w;
		  zgemm_pln.K = n_J;
		  zgemm_pln.N = N_w;

		  zgemm_pln.LDA = N_w;
		  zgemm_pln.LDB = N_w;
		  zgemm_pln.LDC = N_w;

		  zgemm_pln.execute_plan();
		}

		for(int w2=0; w2<N_w; w2++)
		  for(int w1=0; w1<N_w; w1++)
		    M_r_r_w_w(b_i, b_j, r_i, min_r_j, w1, w2) = T_l_times_M_ij_times_T_r[w1 + N_w*w2];
	      }
	    else
	      {
		for(int w2=0; w2<N_w; w2++)
		  for(int w1=0; w1<N_w; w1++)
		    M_r_r_w_w(b_i, b_j, r_i, min_r_j, w1, w2) = 0;
	      }
	  }
	}
      }
    }

    delete [] T;

    delete [] M_ij;
    delete [] T_l;
    delete [] T_r;

    delete [] T_l_times_M_ij;
    delete [] T_l_times_M_ij_times_T_r;

    //     clock_t end = clock();
    //     cout << __FUNCTION__ << "\t" << double(end-begin)/CLOCKS_PER_SEC << endl;
  }

  template<int dimension>
  template<typename configuration_t>
  void cached_nft<dimension>::compute_T(configuration_t& configuration)
  {
//     clock_t begin = clock();

    int N_v = configuration.size();
    int N_w = W.size();

    double x;

    for(int j=0; j<N_v; j++){
      for(int i=0; i<N_w; i++){

	x = (configuration[j].get_tau()*W[i]);

	real(T[i + N_w*j]) = cos(x);
	imag(T[i + N_w*j]) = sin(x);
      }
    }
    
//     clock_t end = clock();
//     cout << __FUNCTION__ << "\t" << double(end-begin)/CLOCKS_PER_SEC << endl;
  }

  template<int dimension>
  template<typename configuration_t>
  int cached_nft<dimension>::sort_configuration(configuration_t&                  configuration,
						std::vector<std::pair<int,int> >& p,
						std::vector<int>&                 start_index,
						std::vector<int>&                 end_index)
  {
//     clock_t begin = clock();

    int N_v = configuration.size();
    static int N_r = DCA_cluster_type::get_cluster_size();

    for(int l=0; l<N_v; l++){
      assert(configuration[l].get_band() == 0);
      p[l].first  = configuration[l].get_r_site();
      p[l].second = l;
    }

//     for(int l=0; l<N_v; l++)
//       cout << p[l].first << "\t" << p[l].second << "\n";
//     cout << "\n";

    sort(p.begin(), p.end());
    
//     for(int l=0; l<N_v; l++)
//       cout << p[l].first << "\t" << p[l].second << "\n";
//     cout << "\n";

    std::vector<int> ordered_r(N_v);
    for(int l=0; l<N_v; l++)
      ordered_r[l] = p[l].first;

    for(int l=0; l<N_r; l++){

      start_index[l] = find(ordered_r.begin(), ordered_r.end(), l) - ordered_r.begin();

      int index = start_index[l];
      while(index < N_v && ordered_r[index] == l)
	index++;
	    
      end_index[l] = index;
    }

//     cout << N_v << endl;
//     for(int l=0; l<N_r; l++)
//       cout << l << "\t" << start_index[l] << "\t" << end_index[l] << "\n";
//     cout << "\n";
    
    int MAX = 0;
    for(int l=0; l<N_r; l++)
      (end_index[l] - start_index[l]) > MAX ? MAX = (end_index[l] - start_index[l]) : MAX = MAX;
     
//     clock_t end = clock();
//     cout << __FUNCTION__ << "\t" << double(end-begin)/CLOCKS_PER_SEC << endl;

    return MAX;
  }

  template<int dimension>
  template<typename configuration_t>
  int cached_nft<dimension>::sort_configuration(configuration_t&                   configuration,
						std::vector<triple<int,int,int> >& p,
						std::vector<int>&                  start_index,
						std::vector<int>&                  end_index)
  {
    static int N_b = b::dmn_size();
    static int N_r = DCA_cluster_type::get_cluster_size();

//     clock_t begin = clock();

    int N_v = configuration.size();
    for(int l=0; l<N_v; l++){
      p[l].first  = configuration[l].get_band();
      p[l].second = configuration[l].get_r_site();
      p[l].third  = l;
    }

    sort(p.begin(),p.end(), triple<int,int,int>::less);

//     for(int l=0; l<N_v; l++)
//       cout << p[l].first << "\t" << p[l].second << "\t" << p[l].third << "\n";
//     cout << "\n";

    for(int b=0; b<N_b; b++){
      for(int r=0; r<N_r; r++){
	triple<int,int,int> t;
	t.first  = b;
	t.second = r;
	
	t.third = 0;
	//cout << lower_bound(p.begin(), p.end(), t, triple<int,int,int>::less)-p.begin() << "\t";
	start_index[b + N_b*r] = lower_bound(p.begin(), p.end(), t, triple<int,int,int>::less)-p.begin();
	  
	t.third = N_v;
	//cout << upper_bound(p.begin(), p.end(), t, triple<int,int,int>::less)-p.begin() << "\n";
	end_index[b + N_b*r] = upper_bound(p.begin(), p.end(), t, triple<int,int,int>::less)-p.begin();
      }
    }

    int MAX = 0;
    for(int l=0; l<N_r; l++)
      (end_index[l] - start_index[l]) > MAX ? MAX = (end_index[l] - start_index[l]) : MAX = MAX;
    
//     clock_t end = clock();
//     cout << __FUNCTION__ << "\t" << double(end-begin)/CLOCKS_PER_SEC << endl;
	
    return MAX;
  }
}


#endif
