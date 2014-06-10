//-*-C++-*-

/** \ingroup DCA */

/*@{*/

/*! \file cached_auxilery_field_values.h  
 *
 *  Contains a templated class over the dimension to represent the cached nfft-plan
 */

#ifndef CACHED_NFFT_H
#define CACHED_NFFT_H

namespace QMC {

  template<int dimension>
  class cached_nfft
  {
  public:
    
    typedef cached_nfft<dimension> this_type;

    const static int    m_accuracy = 6;
    const static double base       = 1.2;
    
    cached_nfft(const cached_nfft<dimension>& obj);
    cached_nfft(int number_of_positive_frequencies);
    ~cached_nfft();

    nfft_plan& get(size_t ceil_log2_k);

    std::complex<double> get_f_hat_at_matsubara_frequency(nfft_plan& nfft, int i);
    std::complex<double> get_f_hat_at_matsubara_frequency(nfft_plan& nfft, int i_0, int i_1);

    int get_number_of_positive_frequencies() const;
    
    void test();
    
    void add_nfft_plans(int number);

  private:

    void initialize_vector();

  private:

    int                    number_of_positive_frequencies;

    size_t                 MAX_VERTICES;

    int                    TOTAL_NFFT_frequencies;
    int*                   NFFT_frequencies;
    int*                   FFTW_frequencies;

    std::vector<nfft_plan> nfft_plan_vector;
  };

  template<int dimension>
  cached_nfft<dimension>::cached_nfft(int number_of_positive_frequencies_input):
    number_of_positive_frequencies(number_of_positive_frequencies_input),
    MAX_VERTICES(1000),
    nfft_plan_vector()
  {
    //cout << __PRETTY_FUNCTION__ << endl;

    TOTAL_NFFT_frequencies = 1;
    NFFT_frequencies = new int[dimension];
    FFTW_frequencies = new int[dimension];

    for(int i=0; i<dimension; i++){
      NFFT_frequencies[i] =   2*2*number_of_positive_frequencies;
      FFTW_frequencies[i] = 2*2*2*number_of_positive_frequencies;
      TOTAL_NFFT_frequencies *= NFFT_frequencies[i];
    }

    initialize_vector();
  }

  template<int dimension>
  cached_nfft<dimension>::cached_nfft(const cached_nfft<dimension>& obj):
    number_of_positive_frequencies(obj.get_number_of_positive_frequencies()),
    MAX_VERTICES(1000),
    nfft_plan_vector()
  {
    //cout << __PRETTY_FUNCTION__ << endl;

    TOTAL_NFFT_frequencies = 1;
    NFFT_frequencies = new int[dimension];
    FFTW_frequencies = new int[dimension];

    for(int i=0; i<dimension; i++){
      NFFT_frequencies[i] =   2*2*number_of_positive_frequencies;
      FFTW_frequencies[i] = 2*2*2*number_of_positive_frequencies;
      TOTAL_NFFT_frequencies *= NFFT_frequencies[i];
    }

    initialize_vector();
  }


  template<int dimension>
  cached_nfft<dimension>::~cached_nfft()
  {
    //cout << __PRETTY_FUNCTION__ << endl;

    for(size_t i=0; i<nfft_plan_vector.size(); i++)
      nfft_finalize(&nfft_plan_vector[i]);

    delete [] NFFT_frequencies;
    delete [] FFTW_frequencies;
  }

  template<int dimension>
  nfft_plan& cached_nfft<dimension>::get(size_t number_of_vertices)
  {
    //cout << "start \t" << __FUNCTION__ << endl;

    int ceil_log_base_k = std::ceil(std::log(double(number_of_vertices))/std::log(base));

    if(ceil_log_base_k < 5)
      ceil_log_base_k = 5;
      
    if(number_of_vertices > MAX_VERTICES)
      add_nfft_plans(number_of_vertices);
     
    // reset
    int nb_nodes = (int)std::pow(base, ceil_log_base_k);

    memset(nfft_plan_vector[ceil_log_base_k].x,     0, nb_nodes*dimension    *sizeof(double) );
    memset(nfft_plan_vector[ceil_log_base_k].f,     0, nb_nodes              *sizeof(std::complex<double>) );
    memset(nfft_plan_vector[ceil_log_base_k].f_hat, 0, TOTAL_NFFT_frequencies*sizeof(std::complex<double>) );

    assert(nb_nodes                == nfft_plan_vector[ceil_log_base_k].M_total);
    assert(TOTAL_NFFT_frequencies  == nfft_plan_vector[ceil_log_base_k].N_total);
    assert(int(number_of_vertices) <= nfft_plan_vector[ceil_log_base_k].M_total);

    //cout << "end \t" << __PRETTY_FUNCTION__ << endl;

    return nfft_plan_vector[ceil_log_base_k];
  }

  template<int dimension>
  void cached_nfft<dimension>::initialize_vector()
  {
    //cout << __FUNCTION__ << endl;

    int MAX_ceil_log_base_k = std::ceil(std::log(double(MAX_VERTICES))/std::log(base))+1;
    nfft_plan_vector.resize(MAX_ceil_log_base_k);

    for(size_t i=0; i<nfft_plan_vector.size(); i++)
      {
	int nb_nodes = (int)std::pow(base,(int)i);

	nfft_init_guru(&nfft_plan_vector[i], dimension, &NFFT_frequencies[0], nb_nodes, &FFTW_frequencies[0], m_accuracy, 
		       MALLOC_X| MALLOC_F_HAT| MALLOC_F| PRE_PHI_HUT| PRE_LIN_PSI| FFTW_INIT| FFT_OUT_OF_PLACE, // nfft-flags
		       FFTW_MEASURE| FFTW_DESTROY_INPUT); // fftw-flags

	nfft_free(nfft_plan_vector[i].psi);
	nfft_plan_vector[i].K= 1024*(m_accuracy+1); // don't even think about touching this !!! 
	nfft_plan_vector[i].psi=(double*) nfft_malloc((nfft_plan_vector[i].K+1)*nfft_plan_vector[i].d*sizeof(double));
	
	nfft_precompute_one_psi(&nfft_plan_vector[i]);
      }

    MAX_VERTICES = (int)std::pow(base,(int)(nfft_plan_vector.size()-1));
  }

  template<int dimension>
  void cached_nfft<dimension>::add_nfft_plans(int number)
  {
    //cout << "start \t" << __FUNCTION__ << endl;

    int MAX_ceil_log_base_k     = std::ceil(std::log(double(number)*base)/std::log(base));
    int old_MAX_ceil_log_base_k = std::ceil(std::log(double(MAX_VERTICES))/std::log(base));

    if(number > int(MAX_VERTICES))
      {
	//cout << number << "\t" << MAX_VERTICES << "\t" << old_MAX_ceil_log_base_k << "\t" <<  MAX_ceil_log_base_k << "\n";

	nfft_plan_vector.resize(MAX_ceil_log_base_k);
	
	for(size_t i=old_MAX_ceil_log_base_k; i<nfft_plan_vector.size(); i++)
	  {	    
	    int nb_nodes = (int)std::pow(base,(int)i);
	    MAX_VERTICES = nb_nodes;
	    
	    nfft_init_guru(&nfft_plan_vector[i], dimension, &NFFT_frequencies[0], nb_nodes, &FFTW_frequencies[0], m_accuracy, 
			   MALLOC_X| MALLOC_F_HAT| MALLOC_F| PRE_PHI_HUT| PRE_LIN_PSI| FFTW_INIT| FFT_OUT_OF_PLACE, // nfft-flags
			   FFTW_MEASURE| FFTW_DESTROY_INPUT); // fftw-flags
	    
// 	    cout << i << "\t";
// 	    cout << nb_nodes << "\t" ;
// 	    cout << m_accuracy << "\t";
// 	    cout << NFFT_frequencies[0] << "\t";
// 	    cout << FFTW_frequencies[0] << "\t";
// 	    cout << nfft_plan_vector[i].M_total << "\t";
// 	    cout << nfft_plan_vector[i].N_total << endl;
	    
	    nfft_free(nfft_plan_vector[i].psi);
	    nfft_plan_vector[i].K= 1024*(m_accuracy+1);
	    nfft_plan_vector[i].psi=(double*) nfft_malloc((nfft_plan_vector[i].K+1)*nfft_plan_vector[i].d*sizeof(double));
	    
	    nfft_precompute_one_psi(&nfft_plan_vector[i]);
	  }

	MAX_VERTICES = (int)std::pow(base,(int)(nfft_plan_vector.size()-1));
      }
  }

  template<int dimension>
  std::complex<double> cached_nfft<dimension>::get_f_hat_at_matsubara_frequency(nfft_plan& nfft, int i)
  {
    throw std::logic_error(__FUNCTION__);
  }

  template<>
  std::complex<double> cached_nfft<1>::get_f_hat_at_matsubara_frequency(nfft_plan& nfft, int i)
  {
    assert(i < NFFT_frequencies[0]/2);
    return std::complex<double>(nfft.f_hat[2*i+1][0], nfft.f_hat[2*i+1][1]);
  }

  template<int dimension>
  std::complex<double> cached_nfft<dimension>::get_f_hat_at_matsubara_frequency(nfft_plan& nfft, int i_0, int i_1)
  {
    throw std::logic_error(__FUNCTION__);
  }

  template<>
  std::complex<double> cached_nfft<2>::get_f_hat_at_matsubara_frequency(nfft_plan& nfft, int i_0, int i_1)
  {
    assert(i_0 < NFFT_frequencies[0]/2); // factor 2 because of even frequencies !!
    assert(i_1 < NFFT_frequencies[1]/2);

    assert((2*i_0+1) < NFFT_frequencies[0]);
    assert((2*i_1+1) < NFFT_frequencies[1]);

    // NFFT is row-major data-alignment !!
    return std::complex<double>(nfft.f_hat[NFFT_frequencies[1]*(2*i_0+1) + (2*i_1+1)][0], 
				nfft.f_hat[NFFT_frequencies[1]*(2*i_0+1) + (2*i_1+1)][1]);

//     return std::complex<double>(nfft.f_hat[(2*i_0+1)+NFFT_frequencies[0]*(2*i_1+1)][0], 
// 				nfft.f_hat[(2*i_0+1)+NFFT_frequencies[0]*(2*i_1+1)][1]);
  }


  template<int dimension>
  int cached_nfft<dimension>::get_number_of_positive_frequencies() const
  {
    return number_of_positive_frequencies;
  }

  template<int dimension>
  void cached_nfft<dimension>::test()
  {
  }

  template<>
  void cached_nfft<1>::test()
  {
    cout << scientific ;
    int size     = 24;
    int nb_nodes = int(std::pow(double(size),1)); 
    
    nfft_plan& nfft = get(nb_nodes);

    for(int i=0; i<size; i++){
	nfft.x[i] = i/double(size)-1./2.;
	nfft.f[i][0]  = sin(2.*3.14* nfft.x[i] );
	nfft.f[i][1]  = 0;

	cout << nfft.f[i][0] << "\t";
      }
      cout << endl;

    nfft_adjoint(&nfft);

    for(int i=0; i<NFFT_frequencies[0]; i++){
 	cout << nfft.f_hat[i] /*<< " ; " << nfft.f_hat[i+NFFT_frequencies[0]*j][1]*/ << "\t";
     }
    cout << endl;

    throw std::logic_error(__FUNCTION__);
  }


  template<>
  void cached_nfft<2>::test()
  {
    cout << scientific ;
    int size     = 24;
    int nb_nodes = int(std::pow(double(size),2)); 
    
    nfft_plan& nfft = get(nb_nodes);

    for(int i=0; i<size; i++){
      for(int j=0; j<size; j++){
	nfft.x[2*(i+size*j)+0] = i/double(size)-1./2.;
	nfft.x[2*(i+size*j)+1] = j/double(size)-1./2.;
	
	nfft.f[i+size*j][0]  = sin(2.*3.1415* (nfft.x[2*(i+size*j)+0] /*+ nfft.x[2*(i+size*j)+1]*/) );
	nfft.f[i+size*j][1]  = 0 ;                // imag

	cout << nfft.f[i+size*j][0] << "\t";
      }
      cout << endl;
    }
    cout << endl;

    nfft_adjoint(&nfft);

    for(int i=0; i<NFFT_frequencies[0]; i++){
      cout << i << "\t";
      for(int j=0; j<NFFT_frequencies[1]; j++){
	cout << nfft.f_hat[i+NFFT_frequencies[0]*j][0]  << "\t";
      }
      cout << endl;
    }
    cout << endl;

    for(int i=0; i<NFFT_frequencies[0]; i++){
      cout << i << "\t";
      for(int j=0; j<NFFT_frequencies[1]; j++){
	cout << nfft.f_hat[i+NFFT_frequencies[0]*j][1]  << "\t";
      }
      cout << endl;
    }
    cout << endl;

    throw std::logic_error(__FUNCTION__);
  }

}


#endif

