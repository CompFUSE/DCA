//-*-C++-*-

#ifndef BASIS_TRANSFORMATIONS_CD_TO_CD_H
#define BASIS_TRANSFORMATIONS_CD_TO_CD_H

namespace MATH_ALGORITHMS
{
  template<typename type_input, typename type_output,int DMN_INDEX>
  class TRANSFORM_DOMAIN<type_input, CONTINUOUS, type_output, CONTINUOUS, DMN_INDEX>
  {
  private:

    const static bool VERBOSE = false;

    typedef basis_transformation<type_input, CONTINUOUS, type_output, CONTINUOUS> contraction_transformation_type;
    typedef typename contraction_transformation_type::matrix_type                 contraction_matrix_type;

    typedef basis_transformation<type_input, CONTINUOUS, type_output, CONTINUOUS> basis_transformation_type;
    typedef typename basis_transformation_type::matrix_type                       matrix_type;

  public:
    
    template<typename scalartype_input, class domain_input, 
	     typename scalartype_output, class domain_output>
    static void execute(FUNC_LIB::function<scalartype_input , domain_input >& f_input, 
			FUNC_LIB::function<scalartype_output, domain_output>& f_output);

    template<typename scalartype, class domain_input, class domain_output>
    static void execute(FUNC_LIB::function<scalartype , domain_input >& f_input, 
			FUNC_LIB::function<scalartype, domain_output>& f_output);

    template<typename scalartype, class domain_input, class domain_output>
    static void execute(FUNC_LIB::function<std::complex<scalartype>, domain_input >& f_input, 
			FUNC_LIB::function<std::complex<scalartype>, domain_output>& f_output);

  private:
 
    template<typename f_input_t, typename f_output_t>
    static void characterize_transformation(f_input_t& f_input, f_output_t f_output, int& M, int& K, int& N, int& P);
    
  };

  template<typename type_input, typename type_output,int DMN_INDEX>
  template<typename scalartype, class domain_input, class domain_output>
  void TRANSFORM_DOMAIN<type_input, CONTINUOUS, type_output, CONTINUOUS, DMN_INDEX>::execute(FUNC_LIB::function<scalartype, domain_input >& f_input, 
											     FUNC_LIB::function<scalartype, domain_output>& f_output)
  {
    
    if(VERBOSE)
      std::cout << "\n\t transform (continuous -> continuous) " << DMN_INDEX << "  " 
	   << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";
    
    int M, K, N, P;
    characterize_transformation(f_input, f_output, M, K, N, P);
    
    matrix_type& T = basis_transformation_type::get_transformation_matrix();
    
    if(VERBOSE)
      {
        std::cout << "\n\t M : " << M << ", K : " << K << ", N : " << N << ", P : " << P << "\n\n"; 
	
	T.print();
	
	f_input.print_fingerprint();
	f_output.print_fingerprint();
      }
    
    for(int l=0; l<P; l++)
      {
	int lin_ind_lhs = M*K*l;
	int lin_ind_rhs = M*N*l;
	
	LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'T', M, N, K,
					     scalartype(1),
					     &f_input(lin_ind_lhs), M,
					     &T(0,0)              , T.get_global_size().first,
					     scalartype(0),
					     &f_output(lin_ind_rhs), M);
      }
  }

  template<typename type_input, typename type_output,int DMN_INDEX>
  template<typename scalartype, class domain_input, class domain_output>
  void TRANSFORM_DOMAIN<type_input, CONTINUOUS, type_output, CONTINUOUS, DMN_INDEX>::execute(FUNC_LIB::function<std::complex<scalartype>, domain_input >& f_input, 
											     FUNC_LIB::function<std::complex<scalartype>, domain_output>& f_output)
  {    
    if(VERBOSE)
      std::cout << "\n\t transform (continuous -> continuous) " << DMN_INDEX << "  " << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";
    
    int M, K, N, P;
    characterize_transformation(f_input, f_output, M, K, N, P);
    
    contraction_matrix_type& T = basis_transformation_type::get_transformation_matrix();
    
    if(VERBOSE)
      {
        std::cout << "\n\t M : " << M << ", K : " << K << ", N : " << N << ", P : " << P << "\n\n"; 
	
	T.print_fingerprint();
	
	f_input.print_fingerprint();
	f_output.print_fingerprint();
      }
    
    scalartype* A = new scalartype[M*K];
    scalartype* C = new scalartype[M*N];

    for(int l=0; l<P; l++)
      {
	int lin_ind_lhs = M*K*l;
	int lin_ind_rhs = M*N*l;

	{// real
	  for(int i=0; i<M*K; i++)
	    A[i] = real(f_input(lin_ind_lhs+i));
	
	  LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'T', M, N, K,
					       scalartype(1),
					       A, M,
					       &T(0,0), T.get_global_size().first,
					       scalartype(0),
					       C, M);

	  for(int i=0; i<M*N; i++)
	    real(f_output(lin_ind_rhs+i)) = C[i];
	}

	{// imag
	  for(int i=0; i<M*K; i++)
	    A[i] = imag(f_input(lin_ind_lhs+i));
	
	  LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'T', M, N, K,
					       scalartype(1),
					       A, M,
					       &T(0,0), T.get_global_size().first,
					       scalartype(0),
					       C, M);

	  for(int i=0; i<M*N; i++)
	    imag(f_output(lin_ind_rhs+i)) = C[i];
	}       	
      }

    delete [] A;
    delete [] C;
  }


  
  template<typename type_input, typename type_output,int DMN_INDEX>
  template<typename f_input_t, typename f_output_t>
  void TRANSFORM_DOMAIN<type_input, CONTINUOUS, type_output, CONTINUOUS, DMN_INDEX>::characterize_transformation(f_input_t& f_input, f_output_t f_output, 
														 int& M, int& K, int& N, int& P)
  {
    M = 1;
    for(int l=0; l<DMN_INDEX; l++)
      M *= f_input[l];
    
    K = f_input [DMN_INDEX];
    N = f_output[DMN_INDEX];
    
    P = 1;
    for(int l=DMN_INDEX+1; l<f_input.signature(); l++)
      P *= f_input[l];
  }
  
}

#endif
