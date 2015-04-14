//-*-C++-*-

#ifndef BASIS_TRANSFORMATIONS_DD_TO_ED_H
#define BASIS_TRANSFORMATIONS_DD_TO_ED_H

#include "fftw3.h"

namespace MATH_ALGORITHMS
{
  template<typename type_input, typename type_output, int DMN_INDEX>
  class TRANSFORM_DOMAIN<type_input, DISCRETE, type_output, EXPANSION, DMN_INDEX> : public TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>
  {
  private:

    const static bool VERBOSE = false;

    typedef typename type_input ::dmn_specifications_type input_specs_type;
    typedef typename type_output::dmn_specifications_type output_specs_type;

    typedef basis_transformation<type_input, DISCRETE, type_output, EXPANSION> basis_transformation_type;
    typedef typename basis_transformation_type::matrix_type                    matrix_type;

  public:

    template<typename scalartype_input, class domain_input, 
	     typename scalartype_output, class domain_output>
    static void execute(FUNC_LIB::function<scalartype_input , domain_input >& f_input, 
			FUNC_LIB::function<scalartype_output, domain_output>& f_output)
    {  
      default_execute(f_input, f_output);
    }

  private:
 
    template<typename scalartype>
    static scalartype vector_norm(scalartype x)
    {
      return abs(x);
    }

    template<typename scalartype>
    static scalartype vector_norm(std::vector<scalartype>& x)
    {
      scalartype result = 0;
      
      for(int l=0; l<x.size(); l++)
	result += x[l]*x[l];

      return result;
    }
    
    static int find_origin()
    {
      int index=0;
      
      for(int l=0; l<type_output::get_size(); l++)
	if(vector_norm(type_output::get_elements()[l])<1.e-6)
	  index = l;
      
      cout << index << "\n";

      return index;
    }

    template<typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
    static void fftw_harmonics_execute(FUNC_LIB::function<scalartype_input , domain_input >& f_input, 
				       FUNC_LIB::function<scalartype_output, domain_output>& f_output,
				       bool renorm);

    template<typename scalartype, class domain_input, class domain_output>
    static void fftw_harmonics_execute(FUNC_LIB::function<std::complex<scalartype>, domain_input >& f_input, 
				       FUNC_LIB::function<std::complex<scalartype>, domain_output>& f_output,
				       bool renorm)
    {
      assert(type_input::dmn_specifications_type::DIMENSION == type_output::dmn_specifications_type::DIMENSION);

      if(VERBOSE)
	cout << "\n\t ifftw-harmonics-transform (discrete -> expansion) " << DMN_INDEX << "  " << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";
      
      int M, K, N, P;
      characterize_transformation(f_input, f_output, M, K, N, P);

      int  rank = type_input::dmn_specifications_type::DIMENSION;
      int* dims = type_input::get_dimensions();
      
      int how_many = M;

      fftw_complex* in  = new fftw_complex[M*K];
      fftw_complex* out = new fftw_complex[M*N];

      int istride = 1;
      int ostride = 1;

      // K=N !
      int idist = K;
      int odist = N; 

      const int* inembed = type_input ::get_dimensions();
      const int* onembed = type_output::get_dimensions();

      if(false and VERBOSE){
	cout << M << "\t" << K << "\t" << N << "\t" << P << "\n";

	cout << rank << "\n";
	for(int i=0; i<rank; i++)
	  cout << dims[i] << "\t";
	cout << "\n";

	f_input .print_fingerprint();
	f_output.print_fingerprint();
      }

      fftw_plan plan = fftw_plan_many_dft(rank, dims, how_many,
					  in , inembed, istride, idist,
					  out, onembed, ostride, odist,
					  FFTW_BACKWARD, FFTW_ESTIMATE);

      int index = find_origin();
      for(int l=0; l<P; l++)      
	{
	  for(int i=0; i<M; i++){     
	    for(int j=0; j<K; j++){     
	      in[j+K*i][0] = real(f_input(M*K*l+i+j*M));
	      in[j+K*i][1] = imag(f_input(M*K*l+i+j*M));
	    }
	  }

	  fftw_execute(plan); 

	  for(int i=0; i<M; i++){     
	    for(int j=0; j<N; j++){     

	      int j_0 = (j-index)<0? (j-index+N) : (j-index);

	      real(f_output(M*N*l+i+j*M)) = out[j_0+N*i][0];
	      imag(f_output(M*N*l+i+j*M)) = out[j_0+N*i][1];
	    }
	  }
	}

      fftw_destroy_plan(plan);

      delete [] in;
      delete [] out;

      if(renorm)
	f_output *= 1./scalartype(K);
    }
  
    template<typename scalartype_input, class domain_input, 
	     typename scalartype_output, class domain_output>
    static void default_execute(FUNC_LIB::function<scalartype_input , domain_input >& f_input, 
				FUNC_LIB::function<scalartype_output, domain_output>& f_output)
    {
      if(VERBOSE)
	cout << "\n default-transform (discrete -> expansion) " << DMN_INDEX << "  " << type_input::get_name() << " --> " << type_output::get_name() << "\n\n";

      matrix_type& T = basis_transformation_type::get_transformation_matrix();

      TRANSFORM_DOMAIN_PROCEDURE<DMN_INDEX>::transform(f_input, f_output, T);
    }
 
  };

}

#endif
