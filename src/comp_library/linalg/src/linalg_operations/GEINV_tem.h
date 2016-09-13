//-*-C++-*-

#ifndef LINALG_GEINV_H
#define LINALG_GEINV_H

namespace LIN_ALG
{

  template<device_type device_t>
  class GEINV
  {
  public:

    template<typename scalartype>
    static void execute(matrix<scalartype, device_t>& A)
    {
      assert(A.size().first  == A.size().second);

      dca::linalg::Vector<int, CPU> IPIV(A.size().first);
      for(int i=0; i<A.size().first; i++)
        IPIV[i] = i+1;

      GETRF<device_t>::execute(A, IPIV.ptr());
      GETRI<device_t>::execute(A, IPIV.ptr());
    }

    template<typename scalartype>
    static inline void execute(int N, scalartype* A)
    {
      int LWORK = 16*std::max(1,N);

      int*        IPIV = new int[N];
      scalartype* WORK = new scalartype[LWORK];
      
      int INFO = -1;
      
      execute(N, A, IPIV, WORK, LWORK, INFO);

      delete [] IPIV;
      delete [] WORK;
    }

    static inline void execute(int N, float* A, int* IPIV, float* WORK, int& LWORK, int& INFO)
    {
      LAPACK::sgetrf_(&N, &N, A, &N, IPIV,               &INFO); assert(INFO==0);
      LAPACK::sgetri_(    &N, A, &N, IPIV, WORK, &LWORK, &INFO); assert(INFO==0);
    }

    static inline void execute(int N, double* A, int* IPIV, double* WORK, int& LWORK, int& INFO)
    {
      LAPACK::dgetrf_(&N, &N, A, &N, IPIV, &INFO);           assert(INFO==0);
      LAPACK::dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO); assert(INFO==0);
    }

    static inline void execute(int N, std::complex<float>* A, int* IPIV, std::complex<float>* WORK, int& LWORK, int& INFO)
    {
      LAPACK::cgetrf_(&N, &N, A, &N, IPIV, &INFO);           assert(INFO==0);
      LAPACK::cgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO); assert(INFO==0);
    }

    static inline void execute(int N, std::complex<double>* A, int* IPIV, std::complex<double>* WORK, int& LWORK, int& INFO)
    {
      LAPACK::zgetrf_(&N, &N, A, &N, IPIV, &INFO);           assert(INFO==0);
      LAPACK::zgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO); assert(INFO==0);
    }

    template<typename scalar_type>
    static bool test(int N, scalar_type* A, scalar_type* A_inv)
    {
      bool OK=true;

      for(int j=0; j<N; ++j){
        for(int i=0; i<N; ++i){

          scalar_type tmp=0;
          for(int l=0; l<N; ++l)
            tmp += A[i+l*N]*A_inv[l+j*N];

          scalar_type res = i==j? 1:0;

          if(abs(tmp-res)>1.e-6)
            OK=false;
        }
      }

      return OK;
    }

    template<typename scalartype>
    static void execute_on_small_matrix(matrix<scalartype, device_t>& A)
    {
      assert(A.size().first == A.size().second);

      switch(A.size().first)
        {
        case 1:
          {
            A(0,0) = 1./A(0,0);
          }
          break;

        default:
          execute(A);
        }
    }

    /*
    template<typename scalartype>
    static void execute_on_Green_function_matrix(matrix<scalartype, device_t>& A)
    {
      assert(A.size().first % 2 == 0);
      assert(A.size().first  == A.size().second);

      switch(A.size().first)
        {
        case 2:
          {
            A(0,0) = scalartype(1.)/A(0,0);
            A(1,1) = scalartype(1.)/A(1,1);
          }
          break;

        default:
          execute(A);
        }
    }
    */
    
    template<typename scalartype>
    class plan
    {
    public:
      
      plan():
	initialized(false),
	N(-1),
	LWORK(-1),
	IPIV(NULL),
	WORK(NULL),
	DATA(NULL)
      {}

      plan(matrix<scalartype, device_t>& A):
	initialized(false),
	N(-1),
	LWORK(-1),
	IPIV(NULL),
	WORK(NULL),
	DATA(NULL)
      {
	initialize(A);
      }

      ~plan()
      {
	if(not (IPIV==NULL)) delete [] IPIV;
	if(not (WORK==NULL)) delete [] WORK;
	if(not (DATA==NULL)) delete [] DATA;
      }
      
      void initialize(matrix<scalartype, device_t>& A)
      {
	assert(A.size().first == A.size().second);
	
	N = A.size().first;

	LWORK = 16*std::max(1,N);

	if(not (IPIV==NULL)) delete [] IPIV;
	if(not (WORK==NULL)) delete [] WORK;
	if(not (DATA==NULL)) delete [] DATA;

	IPIV = new int[N];
	WORK = new scalartype[LWORK];
	DATA = new scalartype[N*N];

	initialized = true;
      }

      void execute(matrix<scalartype, device_t>& A)
      {
	assert(initialized);

	assert(A.size().first == A.size().second);
	assert(A.size().first == N);
	
	for(int j=0; j<N; j++)
	  for(int i=0; i<N; i++)
	    DATA[i+j*N] = A(i,j);

	{
	  int INFO=-1;

	  GEINV::execute(N, DATA, IPIV, WORK, LWORK, INFO);

	  if(INFO!=0) { throw std::logic_error(__FUNCTION__); }
	}

	for(int j=0; j<N; j++)
	  for(int i=0; i<N; i++)
	    A(i,j) = DATA[i+j*N];
      }

    private:
      
      bool initialized;

      int N;

      int LWORK;

      int*        IPIV;
      scalartype* WORK;
      
      scalartype* DATA;
    };

  };
 
}

#endif
