//-*-C++-*-

#ifndef TRANSFORM_DOMAIN_PROCEDURE_H
#define TRANSFORM_DOMAIN_PROCEDURE_H

namespace MATH_ALGORITHMS
{
  template<int DMN_INDEX>
  class TRANSFORM_DOMAIN_PROCEDURE
  {    
    //protected:
  public:

    template<typename f_input_t, typename f_output_t>
    static void characterize_transformation(f_input_t& f_input, f_output_t f_output, int& M, int& K, int& N, int& P);

    template<typename scalartype_1, class domain_input, 
	     typename scalartype_2, class domain_output,
	     typename scalartype_3>
    static void transform(function<scalartype_1, domain_input >&       f_input, 
			  function<scalartype_2, domain_output>&       f_output,
			  LIN_ALG::matrix<scalartype_3, LIN_ALG::CPU>& T);

    template<typename scalartype, class domain_input, class domain_output>
    static void transform(function<scalartype, domain_input >&       f_input, 
			  function<scalartype, domain_output>&       f_output,
			  LIN_ALG::matrix<scalartype, LIN_ALG::CPU>& T);

    template<typename scalartype, class domain_input, class domain_output>
    static void transform(function<std::complex<scalartype>, domain_input >& f_input, 
			  function<std::complex<scalartype>, domain_output>& f_output,
			  LIN_ALG::matrix<scalartype, LIN_ALG::CPU>&         T);

    template<typename scalartype, class domain_input, class domain_output>
    static void transform(function<scalartype, domain_input >&                     f_input, 
			  function<scalartype, domain_output>&                     f_output,
			  LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& T);

    template<class domain_input, class domain_output>
    static void transform(function<float, domain_input >&       f_input, 
			  function<float, domain_output>&       f_output,
			  LIN_ALG::matrix<double, LIN_ALG::CPU>& T);

    template<class domain_input, class domain_output>
    static void transform(function<std::complex<float>, domain_input >&       f_input, 
			  function<std::complex<float>, domain_output>&       f_output,
			  LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU>& T);

    template<typename scalartype, class domain_input, class domain_output>
    static void transform(function<                    scalartype , domain_output>& f_input, 
			  function<       std::complex<scalartype>, domain_input >& f_output,
			  LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU >& T);

    template<typename scalartype, class domain_input, class domain_output>
    static void transform(function<       std::complex<scalartype>, domain_input >& f_input, 
			  function<                    scalartype , domain_output>& f_output,
			  LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU >& T);
    
  };

  template<int DMN_INDEX>
  template<typename f_input_t, typename f_output_t>
  void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::characterize_transformation(f_input_t& f_input, f_output_t f_output, int& M, int& K, int& N, int& P)
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

  template<int DMN_INDEX>
  template<typename scalartype, class domain_input, class domain_output>
  void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(function<scalartype, domain_input >&       f_input, 
							function<scalartype, domain_output>&       f_output,
							LIN_ALG::matrix<scalartype, LIN_ALG::CPU>& T)
  {
    int M, K, N, P;
    characterize_transformation(f_input, f_output, M, K, N, P);

    scalartype alpha(1);
    scalartype beta (0);

    if(M == 1)
      {
	LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'N', T.get_current_size().first, P, T.get_current_size().second,
					     alpha,
					     &T(0,0)     , T.get_global_size().first,
					     &f_input(0) , f_input [DMN_INDEX],
					     beta,
					     &f_output(0), f_output[DMN_INDEX]);
      }
    else
      {
	for(int l=0; l<P; l++)
	  {
	    int lin_ind_lhs = M*K*l;
	    int lin_ind_rhs = M*N*l;
	    
	    LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'T', M, N, K,
						 alpha,
						 &f_input(lin_ind_lhs), M,
						 &T(0,0)              , T.get_global_size().first,
						 beta,
						 &f_output(lin_ind_rhs), M);
	  }
      }
  }

  template<int DMN_INDEX>
  template<typename scalartype, class domain_input, class domain_output>
  void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(function<std::complex<scalartype>, domain_input >& f_input, 
							function<std::complex<scalartype>, domain_output>& f_output,
							LIN_ALG::matrix<scalartype, LIN_ALG::CPU>&         T)
  {
    int M, K, N, P;
    characterize_transformation(f_input, f_output, M, K, N, P);

    for(int l=0; l<P; l++)
      {
	int lin_ind_lhs = M*K*l;
	int lin_ind_rhs = M*N*l;
	
	scalartype alpha(1);
	scalartype beta (0);
	
	LIN_ALG::GEMM<LIN_ALG::CPU>::execute('N', 'T', 2*M, N, K,
					     alpha,
					     &real(f_input(lin_ind_lhs)), 2*M,
					     &T(0,0)              , T.get_global_size().first,
					     beta,
					     &real(f_output(lin_ind_rhs)), 2*M);
      }
  }

  template<int DMN_INDEX>
  template<typename scalartype, class domain_input, class domain_output>
  void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(function<scalartype, domain_input >&                     f_input, 
							function<scalartype, domain_output>&                     f_output,
							LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU>& T)
  {
    LIN_ALG::matrix<scalartype, LIN_ALG::CPU> T_re("T_re", T.get_current_size());
    
    for(int j=0; j<T.get_current_size().second; j++)
      for(int i=0; i<T.get_current_size().first; i++)
	T_re(i,j) = real(T(i,j));

    transform(f_input, f_output, T_re);
  }

  template<int DMN_INDEX>
  template<class domain_input, class domain_output>
  void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(function<       float , domain_input >&       f_input, 
							function<       float , domain_output>&       f_output,
							LIN_ALG::matrix<double, LIN_ALG::CPU>& T)
  {
    LIN_ALG::matrix<float, LIN_ALG::CPU> T_float("T_re", T.get_current_size());

    for(int j=0; j<T.get_current_size().second; j++)
      for(int i=0; i<T.get_current_size().first; i++)
	T_float(i,j) = T(i,j);

    transform(f_input, f_output, T_float);
  }

  template<int DMN_INDEX>
  template<class domain_input, class domain_output>
  void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(function<       std::complex<float> , domain_input >&       f_input, 
							function<       std::complex<float> , domain_output>&       f_output,
							LIN_ALG::matrix<std::complex<double>, LIN_ALG::CPU>& T)
  {
    LIN_ALG::matrix<std::complex<float>, LIN_ALG::CPU> T_float("T_re", T.get_current_size());

    for(int j=0; j<T.get_current_size().second; j++)
      for(int i=0; i<T.get_current_size().first; i++)
	T_float(i,j) = T(i,j);

    transform(f_input, f_output, T_float);
  }

  template<int DMN_INDEX>
  template<typename scalartype, class domain_input, class domain_output>
  void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(function<                    scalartype , domain_output>& f_input, 
							function<       std::complex<scalartype>, domain_input >& f_output,
							LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU >& T)
  {
    function<std::complex<scalartype>, domain_output> f_in("f_in"); 

    for(int i=0; i<f_input.size(); i++)
      real(f_in(i)) = f_input(i);

    transform(f_in, f_output, T);
  }

  template<int DMN_INDEX>
  template<typename scalartype, class domain_input, class domain_output>
  void TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(function<       std::complex<scalartype>, domain_input >& f_input, 
							function<                    scalartype , domain_output>& f_output,
							LIN_ALG::matrix<std::complex<scalartype>, LIN_ALG::CPU >& T)
  {
    f_output = 0.;

    function<scalartype, domain_input > f_in ("f_in");
    function<scalartype, domain_output> f_out("f_out");

    LIN_ALG::matrix<scalartype, LIN_ALG::CPU> T_tmp("T_tmp", T.get_current_size());

    {
      for(int i=0; i<f_input.size(); i++)
	f_in(i) = real(f_input(i));
      
      for(int j=0; j<T.get_current_size().second; j++)
	for(int i=0; i<T.get_current_size().first; i++)
	  T_tmp(i,j) = real(T(i,j));
      
      transform(f_in, f_out, T_tmp);

      for(int i=0; i<f_output.size(); i++)
	f_output(i) += f_out(i);
    }

    {
      for(int i=0; i<f_input.size(); i++)
	f_in(i) = imag(f_input(i));
      
      for(int j=0; j<T.get_current_size().second; j++)
	for(int i=0; i<T.get_current_size().first; i++)
	  T_tmp(i,j) = imag(T(i,j));
      
      transform(f_in, f_out, T_tmp);

      for(int i=0; i<f_output.size(); i++)
	f_output(i) -= f_out(i);
    }

  }

}

#endif
