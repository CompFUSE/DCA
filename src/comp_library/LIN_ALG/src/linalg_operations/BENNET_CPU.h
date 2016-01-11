//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_BENNET_CPU_H
#define LINALG_BENNET_CPU_H

namespace LIN_ALG {

    template<>
    class BENNET<CPU>
    {
    public:

	template<typename scalartype>
	static void execute(matrix<scalartype, CPU>& M, scalartype* c, scalartype* r){
	    
	    assert(M.get_current_size().first == M.get_current_size().second);

	    int N  = M.get_current_size().first;
	    int LD = M.get_global_size() .first;

	    //row_wise_Bennet(N, LD, M.get_ptr(), &c[0], &r[0]);
	    standard_Bennet(N, LD, M.get_ptr(), &c[0], &r[0]);
	}

	template<typename scalartype>
	static void execute_on_row(matrix<scalartype, CPU>& M, int /*c_ind*/, scalartype* /*r*/){
	    
	    assert(M.get_current_size().first == M.get_current_size().second);
            /*
	    int N  = M.get_current_size().first;
	    int LD = M.get_global_size() .first;
	    
	    //row_wise_Bennet(N, LD, M.get_ptr(), &c[0], &r[0]);
	    */
	}

	template<typename scalartype>
	static void execute_on_col(matrix<scalartype, CPU>& M, scalartype* /*c*/, int /*r_ind*/){
	    
	    assert(M.get_current_size().first == M.get_current_size().second);
	    /*
	    int N  = M.get_current_size().first;
	    int LD = M.get_global_size() .first;
	    
	    row_wise_Bennet(N, LD, M.get_ptr(), &c[0], &r[0]);
	    */
	}

	/*!
	  
	  @article{PETER STANGE, ANDREAS GRIEWANK AND MATTHIAS BOLLHO,
	  title = {ON THE EFFICIENT UPDATE OF RECTANGULAR LU FACTORIZATIONS SUBJECT TO LOW RANK MODIFICATIONS},
	  local-url = {file://localhost/Users/peterstaar/Documents/Papers2/Articles/2011/Unknown/2011-11.pdf},
	  }
	*/
	template<typename scalartype>
	static void row_wise_Bennet(int N, int LD, scalartype* M, scalartype* c, scalartype* r){

	    for(int i=0; i<N; i++)
	    {
		for(int j=0; j<i; j++)
		{
		    c[i] = c[i] - c[j]*M[i+j*LD];
		    M[i+j*LD] = M[i+j*LD] + r[j]*c[i];
		}

		M[i+i*LD] = M[i+i*LD] + c[i]*r[i];
		r[i]             = r[i]/M[i+i*LD];
		
		for(int j=i+1; j<N; j++)
		{
		    M[i+j*LD] = M[i+j*LD] + c[i]*r[j];
		    r[j] = r[j] - r[i]*M[i+j*LD];
		}
	    }
	}

	template<typename scalartype>
	static void standard_Bennet(int N, int LD, scalartype* M, scalartype* c, scalartype* r){

	    for(int i=0; i<N; ++i){

		M[i+i*LD] += c[i]*r[i];
		r[i]      /= M[i+i*LD];
		
		for(int j=i+1; j<N; ++j){
		    c[j]      -= c[i]*M[j+i*LD];
		    M[j+i*LD] += c[j]*r[i]; 
		}

		for(int j=i+1; j<N; ++j){
		    M[i+j*LD] += c[i]*r[j]; 
		    r[j]      -= r[i]*M[i+j*LD];
		}
	    }
	}
    };

}

#endif
