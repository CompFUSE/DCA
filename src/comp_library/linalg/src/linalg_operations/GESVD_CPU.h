//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GESVD_CPU_H
#define LINALG_GESVD_CPU_H

namespace LIN_ALG {

  template<>
  class GESVD<CPU>
  {
  public:
    
    /**************************
     ***   general matrix
     **************************/

    template<typename scalartype>
    static void execute(char JOBZ, 
			matrix<scalartype, CPU>& A, 
			vector<scalartype, CPU>& S,
			matrix<scalartype, CPU>& U,
			matrix<scalartype, CPU>& VT);
    
    template<typename scalartype>
    static void execute(char JOBZ, 
			matrix<std::complex<scalartype>, CPU>& A, 
			vector<             scalartype , CPU>& S,
			matrix<std::complex<scalartype>, CPU>& U,
			matrix<std::complex<scalartype>, CPU>& VT);
      
  private:
    
    /**************************
     ***   GESVD-routines
     **************************/
    
    static void execute(char JOBU, char JOBVT, int M, int N, float*  A , int LDA, float*  S, float*  U , int LDU, float*  VT, int LDVT, float*  WORK, int LWORK);
    static void execute(char JOBU, char JOBVT, int M, int N, double* A , int LDA, double* S, double* U , int LDU, double* VT, int LDVT, double* WORK, int LWORK);

    static void execute(char JOBU, char JOBVT, int M, int N, std::complex<float>*  A, int LDA, float*  S, std::complex<float>*  U, int LDU, std::complex<float>*  VT, int LDVT, std::complex<float>*  WORK, int LWORK, float*  RWORK);
    static void execute(char JOBU, char JOBVT, int M, int N, std::complex<double>* A, int LDA, double* S, std::complex<double>* U, int LDU, std::complex<double>* VT, int LDVT, std::complex<double>* WORK, int LWORK, double* RWORK);
    
    /**************************
     ***   GESVDD-routines
     **************************/
    
    static void execute(char JOBZ, int M, int N, float*  A, int LDA, float*  S, float*  U, int LDU, float*  VT, int LDVT, float*  WORK, int LWORK, int* IWORK);
    static void execute(char JOBZ, int M, int N, double* A, int LDA, double* S, double* U, int LDU, double* VT, int LDVT, double* WORK, int LWORK, int* IWORK);
    
    static void execute(char JOBZ, int M, int N, std::complex<float>*  A, int LDA, float*  S, std::complex<float>*  U, int LDU, std::complex<float>*  VT, int LDVT, std::complex<float>*  WORK, int LWORK, float*  RWORK, int* IWORK);
    static void execute(char JOBZ, int M, int N, std::complex<double>* A, int LDA, double* S, std::complex<double>* U, int LDU, std::complex<double>* VT, int LDVT, std::complex<double>* WORK, int LWORK, double* RWORK, int* IWORK);
  };
  
  template<typename scalartype>
  void GESVD<CPU>::execute(char JOBZ, 
			   matrix<scalartype, CPU>& A, 
			   vector<scalartype, CPU>& S,
			   matrix<scalartype, CPU>& U,
			   matrix<scalartype, CPU>& VT)
  {
    if( JOBZ != 'N' and JOBZ != 'A')
      throw std::logic_error(__FUNCTION__);
    
    matrix<scalartype, CPU> X;
    X.copy_from(A);
    
    int M = X.get_current_size().first;
    int N = X.get_current_size().second;
    
    int LDX = X.get_global_size().first;
    
    int LDU  = U .get_global_size().first;
    int LDVT = VT.get_global_size().first;

    if(false)
      {
	char JOBU  = JOBZ;
	char JOBVT = JOBZ;

	int LWORK  = -1;

	vector<scalartype, CPU> WORK(1);

	{// optimal work-space query
	  execute(JOBU, JOBVT, M, N, &X(0,0), LDX, &S[0], &U(0,0), LDU, &VT(0,0), LDVT, &WORK[0], LWORK);

	  LWORK = WORK[0];
	  WORK.resize(LWORK);
	}

	execute(JOBU, JOBVT, M, N, &X(0,0), LDX, &S[0], &U(0,0), LDU, &VT(0,0), LDVT, &WORK[0], LWORK);
      }
    else
      {
	int LWORK  = -1;
		
	int LIWORK = 8*std::min(M,N);
	
	vector<scalartype, CPU> WORK(1);
	vector<int       , CPU> IWORK(LIWORK);
	
	{// optimal work-space query
	  execute(JOBZ, M, N, &X(0,0), LDX, &S[0], &U(0,0), LDU, &VT(0,0), LDVT, &WORK[0], LWORK, &IWORK[0]);
	  
	  LWORK = WORK[0];
	  WORK.resize(LWORK);
	}

	execute(JOBZ, M, N, &X(0,0), LDX, &S[0], &U(0,0), LDU, &VT(0,0), LDVT, &WORK[0], LWORK, &IWORK[0]);	  
      }
  }
  
  template<typename scalartype>
  void GESVD<CPU>::execute(char JOBZ, 
			   matrix<std::complex<scalartype>, CPU>& A, 
			   vector<             scalartype , CPU>& S,
			   matrix<std::complex<scalartype>, CPU>& U,
			   matrix<std::complex<scalartype>, CPU>& VT)
  {
    if( JOBZ != 'N' and JOBZ != 'A')
      throw std::logic_error(__FUNCTION__);
    
    matrix<std::complex<scalartype>, CPU> X;
    X.copy_from(A);
    
    int M = X.get_current_size().first;
    int N = X.get_current_size().second;
    
    int LDX = X.get_global_size().first;
    
    int LDU  = U .get_global_size().first;
    int LDVT = VT.get_global_size().first;

    if(false)
      {
	char JOBU  = JOBZ;
	char JOBVT = JOBZ;

	int LWORK  = -1;
	int LRWORK = 5*std::min(M,N);

	vector<std::complex<scalartype>, CPU> WORK(1);
	vector<             scalartype , CPU> RWORK(LRWORK);

	{// optimal work-space query
	  execute(JOBU, JOBVT, M, N, &X(0,0), LDX, &S[0], &U(0,0), LDU, &VT(0,0), LDVT, &WORK[0], LWORK, &RWORK[0]);
	  
	  LWORK = real(WORK[0]);
	  WORK.resize(LWORK);
	}

	execute(JOBU, JOBVT, M, N, &X(0,0), LDX, &S[0], &U(0,0), LDU, &VT(0,0), LDVT, &WORK[0], LWORK, &RWORK[0]);
      }
    else
      {
	int LWORK  = -1;
	int LRWORK = -1;
	
	if(JOBZ == 'N') 
	  LRWORK = 7*std::min(M,N);
	else
	  LRWORK = 5*std::min(M,N)*std::min(M,N) + 5*std::min(M,N);
	
	int LIWORK = 8*std::min(M,N);
	
	vector<std::complex<scalartype>, CPU> WORK(1);
	vector<             scalartype , CPU> RWORK(LRWORK);
	vector<             int        , CPU> IWORK(LIWORK);
	
	{// optimal work-space query
	  execute(JOBZ, M, N, &X(0,0), LDX, &S[0], &U(0,0), LDU, &VT(0,0), LDVT, &WORK[0], LWORK, &RWORK[0], &IWORK[0]);
	  
	  //cout << "\n\t optimal work-space for SVD : " << real(WORK[0]) << "\n";

	  LWORK = real(WORK[0]);
	  WORK.resize(LWORK);
	}

	execute(JOBZ, M, N, &X(0,0), LDX, &S[0], &U(0,0), LDU, &VT(0,0), LDVT, &WORK[0], LWORK, &RWORK[0], &IWORK[0]);	  
      }
  }
  

  /**************************
   ***   GESVD-routines
   **************************/

  void GESVD<CPU>::execute(char JOBU, char JOBVT, int M, int N, float*  A , int LDA, float*  S, float*  U , int LDU, float*  VT, int LDVT, float*  WORK, int LWORK)
  {
    int INFO = -1;

    LAPACK::sgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);
    
    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }
  
  void GESVD<CPU>::execute(char JOBU, char JOBVT, int M, int N, double*  A , int LDA, double*  S, double*  U , int LDU, double*  VT, int LDVT, double*  WORK, int LWORK)
  {
    int INFO = -1;

    LAPACK::dgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);
    
    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  void GESVD<CPU>::execute(char JOBU, char JOBVT, int M, int N, std::complex<float>*  A, int LDA, float*  S, std::complex<float>*  U, int LDU, std::complex<float>*  VT, int LDVT, std::complex<float>*  WORK, int LWORK, float*  RWORK)
  {
    int INFO = -1;

    LAPACK::cgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, RWORK, &INFO);
    
    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }
  
  void GESVD<CPU>::execute(char JOBU, char JOBVT, int M, int N, std::complex<double>* A, int LDA, double* S, std::complex<double>* U, int LDU, std::complex<double>* VT, int LDVT, std::complex<double>* WORK, int LWORK, double* RWORK)
  {
    int INFO = -1;

    LAPACK::zgesvd_(&JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, RWORK, &INFO);
    
    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  /**************************
   ***   GESVDD-routines
   **************************/

  void GESVD<CPU>::execute(char JOBZ, int M, int N, float*  A, int LDA, float*  S, float*  U, int LDU, float*  VT, int LDVT, float* WORK, int LWORK, int* IWORK)
  {
    int INFO = -1;

    LAPACK::sgesdd_(&JOBZ,&M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &INFO);
    
    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }
  
  void GESVD<CPU>::execute(char JOBZ, int M, int N, double* A, int LDA, double* S, double* U, int LDU, double* VT, int LDVT, double* WORK, int LWORK, int* IWORK)
  {
    int INFO = -1;

    LAPACK::dgesdd_(&JOBZ,&M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &INFO);
    
    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }
  
  void GESVD<CPU>::execute(char JOBZ, int M, int N, std::complex<float>*  A, int LDA, float*  S, std::complex<float>*  U, int LDU, std::complex<float>*  VT, int LDVT, std::complex<float>*  WORK, int LWORK, float*  RWORK, int* IWORK)
  {
    int INFO = -1;

    LAPACK::cgesdd_(&JOBZ,&M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, RWORK, IWORK, &INFO);
    
    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }
  
  void GESVD<CPU>::execute(char JOBZ, int M, int N, std::complex<double>* A, int LDA, double* S, std::complex<double>* U, int LDU, std::complex<double>* VT, int LDVT, std::complex<double>* WORK, int LWORK, double* RWORK, int* IWORK)
  {
    int INFO = -1;

    LAPACK::zgesdd_(&JOBZ,&M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, RWORK, IWORK, &INFO);
    
    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

}

#endif

