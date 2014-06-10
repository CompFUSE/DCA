//-*-C++-*-

/** \ingroup DCA */

/*@{*/

/*! \file cached_auxilery_field_values.h  
 *
 *  Contains a templated class over the dimension to represent the cached nnfft-plan
 */

#ifndef CACHED_NNFFT_H
#define CACHED_NNFFT_H

namespace QMC {

  template<int dimension>
  class cached_nnfft
  {
#include "type_definitions.h"

  public:
    
    typedef cached_nnfft<dimension> this_type;

    const static int    m_accuracy = 4;
    const static double base       = 1.2;
    
    cached_nnfft();
    cached_nnfft(const cached_nnfft<dimension>& obj);
    ~cached_nnfft();

    nnfft_plan& get(size_t ceil_log2_k);

    std::complex<double> get_f_hat_at_matsubara_frequency(nnfft_plan& nnfft, int i);
    std::complex<double> get_f_hat_at_matsubara_frequency(nnfft_plan& nnfft, int i_0, int i_1);

    int get_number_of_positive_frequencies() const;
    
    static void sinus_test();
    static void test();
    
    void add_nnfft_plans(int number);

  private:

    void initialize_vector();

  private:

    std::vector<double> matsubara_frequencies;

    int                 number_of_positive_frequencies;
    int                 number_of_NNFFT_positive_frequencies;

    size_t              MAX_VERTICES;

    int                 N_total;
    int*                NNFFT_frequencies;
    int*                FFTW_frequencies;
    
    std::vector<nnfft_plan> nnfft_plan_vector;
  };

  template<int dimension>
  cached_nnfft<dimension>::cached_nnfft():
    matsubara_frequencies               (frequency_vertex_domain_type::get_elements()),
    number_of_positive_frequencies      (    matsubara_frequencies.size()/2),
    number_of_NNFFT_positive_frequencies(int(matsubara_frequencies.back()*DCA_parameters_type::beta/M_PI)+1),
    MAX_VERTICES(1000),
    nnfft_plan_vector()
  {
    N_total = 1;

    NNFFT_frequencies = new int[dimension];
    FFTW_frequencies  = new int[dimension];

    for(int i=0; i<dimension; i++){
      NNFFT_frequencies[i] =   2*number_of_NNFFT_positive_frequencies;
      FFTW_frequencies[i]  = 2*2*number_of_NNFFT_positive_frequencies;
      N_total             *= (2*number_of_positive_frequencies);
    }

    initialize_vector();
  }

  template<int dimension>
  cached_nnfft<dimension>::cached_nnfft(const cached_nnfft<dimension>& obj):
    matsubara_frequencies               (frequency_vertex_domain_type::get_elements()),
    number_of_positive_frequencies      (    matsubara_frequencies.size()/2),
    number_of_NNFFT_positive_frequencies(int(matsubara_frequencies.back()*DCA_parameters_type::beta/M_PI)),
    MAX_VERTICES(1000),
    nnfft_plan_vector()
  {
    N_total = 1;
    NNFFT_frequencies = new int[dimension];
    FFTW_frequencies = new int[dimension];

    for(int i=0; i<dimension; i++){
      NNFFT_frequencies[i] =   2*number_of_positive_frequencies;
      FFTW_frequencies[i]  = 2*2*number_of_positive_frequencies;
      N_total             *= NNFFT_frequencies[i];
    }

    initialize_vector();
  }


  template<int dimension>
  cached_nnfft<dimension>::~cached_nnfft()
  {
    //cout << __PRETTY_FUNCTION__ << endl;

    for(size_t i=0; i<nnfft_plan_vector.size(); i++)
      nnfft_finalize(&nnfft_plan_vector[i]);

    delete [] NNFFT_frequencies;
    delete [] FFTW_frequencies;
  }

  template<int dimension>
  nnfft_plan& cached_nnfft<dimension>::get(size_t number_of_vertices)
  {
    //cout << "start \t" << __FUNCTION__ << endl;

    if(number_of_vertices > MAX_VERTICES)
      add_nnfft_plans(number_of_vertices);

    int ceil_log_base_k = std::ceil(std::log(double(number_of_vertices))/std::log(base));

    if(ceil_log_base_k < 6)
      ceil_log_base_k = 6;
           
    // reset
    int nb_nodes = nnfft_plan_vector[ceil_log_base_k].M_total;
    int nb_freq  = nnfft_plan_vector[ceil_log_base_k].N_total;

    memset(nnfft_plan_vector[ceil_log_base_k].x,     0, nb_nodes*dimension*sizeof(double) );
    memset(nnfft_plan_vector[ceil_log_base_k].f,     0, nb_nodes*sizeof(std::complex<double>) );
    memset(nnfft_plan_vector[ceil_log_base_k].f_hat, 0, nb_freq *sizeof(std::complex<double>) );

    assert(int(number_of_vertices) <= nnfft_plan_vector[ceil_log_base_k].M_total);

    return nnfft_plan_vector[ceil_log_base_k];
  }


  template<int dimension>
  void cached_nnfft<dimension>::initialize_vector()
  {
    //cout << __FUNCTION__ << endl;

    int MAX_ceil_log_base_k = std::ceil(std::log(double(MAX_VERTICES))/std::log(base))+1;
    nnfft_plan_vector.resize(MAX_ceil_log_base_k);

    for(size_t i=0; i<nnfft_plan_vector.size(); i++)
      {
	int M_total = int(ceil(std::pow(base, int(i))));

	nnfft_init_guru(&nnfft_plan_vector[i], dimension, N_total, M_total, NNFFT_frequencies, FFTW_frequencies, m_accuracy, 
			PRE_LIN_PSI| PRE_PHI_HUT| MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F); // nfft-flags
	
	/** init (uneven) freq-nodes */
	switch(dimension)
	{
	case 1:
	  {
	    for(int l=0; l<N_total; l++)
	      for(int d=0;d<dimension;d++)
		nnfft_plan_vector[i].v[d + dimension*l] = matsubara_frequencies[l]*(DCA_parameters_type::beta/M_PI)/double(NNFFT_frequencies[0]);
	  }
	  break;

	case 2:
	  {
	    assert(size_t(sqrt(N_total)) == matsubara_frequencies.size());
	    
	    int N_d = matsubara_frequencies.size();
	    
	    for(int l1=0; l1<N_d; l1++){
	      for(int l2=0; l2<N_d; l2++){
		nnfft_plan_vector[i].v[0 + dimension*(l1+N_d*l2)] = matsubara_frequencies[l1]*(DCA_parameters_type::beta/M_PI)/double(NNFFT_frequencies[0]);
		nnfft_plan_vector[i].v[1 + dimension*(l1+N_d*l2)] = matsubara_frequencies[l2]*(DCA_parameters_type::beta/M_PI)/double(NNFFT_frequencies[1]);
	      }
	    }
	  }
	  break;

	default:
	  throw std::logic_error(__FUNCTION__);
	}

	/** precompute psi, the entries of the matrix B */
	nnfft_precompute_lin_psi(&nnfft_plan_vector[i]);

	MAX_VERTICES = M_total;
      }
  }

  template<int dimension>
  void cached_nnfft<dimension>::add_nnfft_plans(int number)
  {
    //cout << "start \t" << __FUNCTION__ << endl;

    int MAX_ceil_log_base_k     = std::ceil(std::log(double(number)*base)/std::log(base));

    int old_size = nnfft_plan_vector.size();

    if(number > int(MAX_VERTICES))
      {
	nnfft_plan_vector.resize(MAX_ceil_log_base_k);
	
	for(size_t i=old_size; i<nnfft_plan_vector.size(); i++)
	  {	    
	    int M_total = int(ceil(std::pow(base, int(i))));

	    nnfft_init_guru(&nnfft_plan_vector[i], dimension, N_total, M_total, NNFFT_frequencies, FFTW_frequencies, m_accuracy, 
			    PRE_LIN_PSI| PRE_PHI_HUT| MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F); // nfft-flags
	
	    /** init (uneven) freq-nodes */
	    switch(dimension)
	      {
	      case 1:
		{
		  for(int l=0; l<N_total; l++)
		    for(int d=0;d<dimension;d++)
		      nnfft_plan_vector[i].v[d + dimension*l] = matsubara_frequencies[l]*(DCA_parameters_type::beta/M_PI)/double(NNFFT_frequencies[0]);
		}
		break;

	      case 2:
		{
		  assert(size_t(sqrt(N_total)) == matsubara_frequencies.size());
	    
		  int N_d = matsubara_frequencies.size();
	    
		  for(int l1=0; l1<N_d; l1++){
		    for(int l2=0; l2<N_d; l2++){
		      nnfft_plan_vector[i].v[0 + dimension*(l1+N_d*l2)] = matsubara_frequencies[l1]*(DCA_parameters_type::beta/M_PI)/double(NNFFT_frequencies[0]);
		      nnfft_plan_vector[i].v[1 + dimension*(l1+N_d*l2)] = matsubara_frequencies[l2]*(DCA_parameters_type::beta/M_PI)/double(NNFFT_frequencies[1]);
		    }
		  }
		}
		break;

	      default:
		throw std::logic_error(__FUNCTION__);
	      }

	    /** precompute psi, the entries of the matrix B */
	    nnfft_precompute_lin_psi(&nnfft_plan_vector[i]);

	    MAX_VERTICES = M_total;
	  }	
      }

    assert(number <= int(MAX_VERTICES));
  }

  template<int dimension>
  std::complex<double> cached_nnfft<dimension>::get_f_hat_at_matsubara_frequency(nnfft_plan& nnfft, int i)
  {
    throw std::logic_error(__FUNCTION__);
  }

  template<>
  std::complex<double> cached_nnfft<1>::get_f_hat_at_matsubara_frequency(nnfft_plan& nnfft, int i)
  {
    assert(i >= 0 && size_t(i) < matsubara_frequencies.size());
    return std::complex<double>(nnfft.f_hat[i][0], nnfft.f_hat[i][1]);
  }

  template<int dimension>
  std::complex<double> cached_nnfft<dimension>::get_f_hat_at_matsubara_frequency(nnfft_plan& nnfft, int i_0, int i_1)
  {
    throw std::logic_error(__FUNCTION__);
  }

  template<>
  std::complex<double> cached_nnfft<2>::get_f_hat_at_matsubara_frequency(nnfft_plan& nnfft, int i_0, int i_1)
  {
    assert(i_0 >= 0 && size_t(i_0) < matsubara_frequencies.size());
    assert(i_1 >= 0 && size_t(i_1) < matsubara_frequencies.size());

    return std::complex<double>(nnfft.f_hat[i_0 + matsubara_frequencies.size()*i_1][0], 
				nnfft.f_hat[i_0 + matsubara_frequencies.size()*i_1][1]);

  }

  template<int dimension>
  void cached_nnfft<dimension>::test()
  {
    cout << __PRETTY_FUNCTION__ << endl;
  }

  template<>
  void cached_nnfft<1>::test()
  {
    cout << __PRETTY_FUNCTION__ << endl;
  }

  template<>
  void cached_nnfft<2>::test()
  {
    cout << __PRETTY_FUNCTION__ << endl;

    cout << scientific ;

    double a1 = 5;
    double a2 = 1;

    //size_t N_d = frequency_vertex_domain_type::get_elements().size();

    QMC::cached_nnfft<2> nnfft_obj;

    for(size_t M=2; M<6; M+=1)
      {
	int M_total = M*M;

	nnfft_plan& nnfft_pln = nnfft_obj.get(M_total);

	cout << nnfft_pln.M_total << "\t" << M_total << "\t";

	clock_t t = clock();

	for(size_t j1=0;j1<M;j1++)
	  {
	    for(size_t j2=0;j2<M;j2++)
	      {
		nnfft_pln.x[0 + 2*(j1 + M*j2)]=double(j1)/double(M)-0.5; 
		nnfft_pln.x[1 + 2*(j1 + M*j2)]=double(j2)/double(M)-0.5;
	      }
	  }
	cout << double(clock()-t)/double(CLOCKS_PER_SEC) << "\t";

// 	for(size_t j1=0;j1<M;j1++)
// 	  for(size_t j2=0;j2<M;j2++)
// 	    cout << nnfft_pln.x[0 + 2*(j1 + M*j2)] << "\t" << nnfft_pln.x[1 + 2*(j1 + M*j2)] << endl;

// 	for(size_t j1=0;j1<N_d;j1++)
// 	  for(size_t j2=0;j2<N_d;j2++)
// 	    cout << nnfft_pln.v[0 + 2*(j1 + N_d*j2)] << "\t" << nnfft_pln.v[1 + 2*(j1 + N_d*j2)] << endl;


	//cout << "precompute NNFFT \t" ;
	clock_t t1 = clock();
	nnfft_precompute_phi_hut(&nnfft_pln);
	cout << double(clock()-t1)/double(CLOCKS_PER_SEC) << "\t";

	for(size_t j1=0;j1<M;j1++)
	  {
	    for(size_t j2=0;j2<M;j2++)
	      {
		nnfft_pln.f[j1 + M*j2][0] = sin(2*M_PI*(a1*nnfft_pln.x[0+2*(j1 + M*j2)] + a2*nnfft_pln.x[1+2*(j1 + M*j2)]));
		nnfft_pln.f[j1 + M*j2][1] = 0;
	      }
	  }

// 	for(size_t l1=0; l1<M; l1++){
// 	  for(size_t l2=0; l2<M; l2++){
// 	    cout << nnfft_pln.f[(l1 + M*l2)][0] << "\t";
// 	  }
// 	  cout << "\n";
// 	}
// 	cout << "\n";

	//cout << "NNFFT \t" ;
	clock_t t2 = clock();
	nnfft_adjoint(&nnfft_pln);
	cout << double(clock()-t2)/double(CLOCKS_PER_SEC) << "\t";
	
// 	for(size_t l1=0; l1<N_d; l1++){
// 	  for(size_t l2=0; l2<N_d; l2++){
// 	    cout << nnfft_pln.f_hat[l1 + N_d*l2][1] << "\t";
// 	  }
// 	  cout << "\n";
// 	}
// 	cout << "\n";
	
	cout << "\n";
      }

    throw std::logic_error(__FUNCTION__);
  }

  template<int dimension>
  void cached_nnfft<dimension>::sinus_test()
  {
  }

  template<>
  void cached_nnfft<1>::sinus_test()
  {
    cout << scientific ;

    int a = 2;

    int j;                              /**< index for nodes and freqencies   */
    nnfft_plan my_plan;                 /**< plan for the nfft                */
    
    int N[1];
    N[0]=20;

    /** init an one dimensional plan */
    nnfft_init(&my_plan, 1, 20, 40, N);
    
    /** init pseudo random nodes */
    for(j=0;j<my_plan.M_total;j++)
      {
	my_plan.x[j]=double(j)/double(my_plan.M_total)-0.5;
      }
    /** init pseudo random nodes */
    for(j=0;j<my_plan.N_total;j++)
      {
	my_plan.v[j]=double(j)/double(N[0])-double(my_plan.N_total/2)/double(N[0]);
      }
    
    /** precompute psi, the entries of the matrix B */
    if(my_plan.nnfft_flags & PRE_PSI)
      nnfft_precompute_psi(&my_plan);
    
    if(my_plan.nnfft_flags & PRE_FULL_PSI)
      nnfft_precompute_full_psi(&my_plan);

    if(my_plan.nnfft_flags & PRE_LIN_PSI)
      nnfft_precompute_lin_psi(&my_plan);
    
    /** precompute phi_hut, the entries of the matrix D */
    if(my_plan.nnfft_flags & PRE_PHI_HUT)
      nnfft_precompute_phi_hut(&my_plan);
    
    /** init pseudo random Fourier coefficients and show them */
    for(j=0;j<my_plan.M_total;j++){
      my_plan.f[j][0] = sin(2*M_PI*a*my_plan.x[j]);
      my_plan.f[j][1] = 0;

      cout << my_plan.x[j] << "\t" << my_plan.f[j][0] << "\t" << my_plan.f[j][1] << endl;
    }
    cout << endl;
    
    /** direct trafo and show the result */
    nndft_adjoint(&my_plan);
    for(j=0;j<my_plan.N_total;j++)
      cout << my_plan.v[j] << /*"\t" << my_plan.f_hat[j][0] <<*/ "\t" << my_plan.f_hat[j][1] << endl;
     cout << endl;
    
    /** approx. trafo and show the result */
    nnfft_adjoint(&my_plan);
    for(j=0;j<my_plan.N_total;j++)
      cout << my_plan.v[j] << /*"\t" << my_plan.f_hat[j][0] <<*/ "\t" << my_plan.f_hat[j][1] << endl;
    cout << endl;
    
    /** finalise the one dimensional plan */
    nnfft_finalize(&my_plan);


    throw std::logic_error(__FUNCTION__);
  }


  template<>
  void cached_nnfft<2>::sinus_test()
  {
    cout << __PRETTY_FUNCTION__ << endl;

    cout << scientific ;

    int a1 = 5;
    int a2 = 1;

    //    int j;                              /**< index for nodes and freqencies   */
    nnfft_plan my_plan;                 /**< plan for the nfft                */
    
    int dimension = 2;

    int M[dimension];
    M[0] = 40;
    M[1] = 40;

    int M_total = 1.;
    for(int l=0; l<dimension; l++)
      M_total *= M[l];

    int N_w[dimension];
    N_w[0] = 20;
    N_w[1] = 20;

    int N_total = 1.;
    for(int l=0; l<dimension; l++)
      N_total *= N_w[l];

    int NFFT_frequencies[dimension];
    NFFT_frequencies[0]= 100; // max positive frequency
    NFFT_frequencies[1]= 100; // max positive frequency

    int FFTW_frequencies[dimension];
    FFTW_frequencies[0]= 2*NFFT_frequencies[0]; // max positive frequency
    FFTW_frequencies[1]= 2*NFFT_frequencies[1]; // max positive frequency

    /** init an one dimensional plan */
    //     nnfft_init(&my_plan, 2, N_total, M_total, NFFT_frequencies);
    nnfft_init_guru(&my_plan, dimension, N_total, M_total, NFFT_frequencies, FFTW_frequencies, m_accuracy, 
		    PRE_LIN_PSI| PRE_PHI_HUT| MALLOC_X| MALLOC_V| MALLOC_F_HAT| MALLOC_F); // nfft-flags
    
    /** init (uneven) freq-nodes */
    {
      for(int j1=0;j1<N_w[0];j1++)
	{
	  for(int j2=0;j2<N_w[1];j2++)
	    {
	      my_plan.v[0 + 2*(j1 + N_w[0]*j2)] = double(2*(j1-N_w[0]/2)+1)/double(NFFT_frequencies[0]);
	      my_plan.v[1 + 2*(j1 + N_w[0]*j2)] = double(2*(j2-N_w[1]/2)+1)/double(NFFT_frequencies[1]);
	    }
	}
    }

    /** precompute psi, the entries of the matrix B */
    clock_t t2 = clock();
    if(my_plan.nnfft_flags & PRE_LIN_PSI)
      nnfft_precompute_lin_psi(&my_plan);
    cout << double(clock()-t2)/double(CLOCKS_PER_SEC) << endl;

    /** init nodes */
    for(int j1=0;j1<M[0];j1++)
      {
	for(int j2=0;j2<M[1];j2++)
	  {
	    my_plan.x[0 + 2*(j1 + M[0]*j2)]=double(j1)/double(M[0])-0.5;
	    my_plan.x[1 + 2*(j1 + M[0]*j2)]=double(j2)/double(M[1])-0.5;
	  }
      }

    cout << "start precompute ... " << endl;

    /** precompute phi_hut, the entries of the matrix D */
    clock_t t1 = clock();
    if(my_plan.nnfft_flags & PRE_PHI_HUT)
      nnfft_precompute_phi_hut(&my_plan);
    cout << double(clock()-t1)/double(CLOCKS_PER_SEC) << endl;
    
    /** precompute psi, the entries of the matrix B */
//     if(my_plan.nnfft_flags & PRE_PSI)
//       nnfft_precompute_psi(&my_plan);
    
//     if(my_plan.nnfft_flags & PRE_FULL_PSI)
//       nnfft_precompute_full_psi(&my_plan);

//     if(my_plan.nnfft_flags & PRE_LIN_PSI)
//       nnfft_precompute_lin_psi(&my_plan);
    
//     /** precompute phi_hut, the entries of the matrix D */
//     if(my_plan.nnfft_flags & PRE_PHI_HUT)
//       nnfft_precompute_phi_hut(&my_plan);
    
    
    /** init pseudo random Fourier coefficients and show them */
    for(int j1=0;j1<M[0];j1++)
      {
	for(int j2=0;j2<M[1];j2++)
	  {
	    my_plan.f[j1 + M[0]*j2][0] = sin(2*M_PI*(a1*my_plan.x[0+2*(j1 + M[0]*j2)] + a2*my_plan.x[1+2*(j1 + M[0]*j2)]));
	    my_plan.f[j1 + M[0]*j2][1] = 0;
	    
	    cout /*<< my_plan.x[0+2*(j1 + M[0]*j2)] << "\t" << my_plan.x[1+2*(j1 + M[0]*j2)]*/
	      << "\t" 
	      << my_plan.f[j1 + M[0]*j2][0] /*<< "\t" << my_plan.f[j1 + M[0]*j2][1] << "\n"*/;
	  }
	cout << endl;
      }
    cout << endl;
    
    /** direct trafo and show the result */
    cout << "NNDFT" << endl;
    nndft_adjoint(&my_plan);

    /** approx. trafo and show the result */
    cout << "NNFFT \t" ;
    clock_t t = clock();
    nnfft_adjoint(&my_plan);
    cout << double(clock()-t)/double(CLOCKS_PER_SEC) << endl;

    for(int j1=0;j1<N_w[0];j1++){
      for(int j2=0;j2<N_w[1];j2++)
	cout /*<< my_plan.v[0 + 2*(j1 + M[0]*j2)] << "\t" << my_plan.v[1 + 2*(j1 + N_w[0]*j2)]*/ << "\t" << my_plan.f_hat[j1 + N_w[0]*j2][1];
      cout << endl;
    }
    cout << endl;

    /** finalise the one dimensional plan */
    nnfft_finalize(&my_plan);

    throw std::logic_error(__FUNCTION__);
  }


}

#endif

