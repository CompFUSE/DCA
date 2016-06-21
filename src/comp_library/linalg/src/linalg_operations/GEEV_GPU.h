 //-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

 #ifndef LINALG_GEEV_GPU_H
 #define LINALG_GEEV_GPU_H

 namespace LIN_ALG {

   namespace GPU_KERNELS_GEEV
   {
     /**************************
      ***   GEEV-routines
      **************************/

     void sgeev(char JOBVL,  char JOBVR, int N, float* A,  int LDA, float* WR, float* WI, float* VL, int LDVL, float* VR, int LDVR, float* WORK, int LWORK);
     void dgeev(char JOBVL,  char JOBVR, int N, double* A,  int LDA, double* WR, double* WI, double* VL, int LDVL, double* VR, int LDVR, double* WORK, int LWORK);
     void cgeev(char JOBVL,  char JOBVR, int N, std::complex<float>* A,  int LDA, std::complex<float>* W, std::complex<float>* VL, int LDVL, std::complex<float>* VR, int LDVR, std::complex<float>* WORK, int LWORK, float* RWORK);
     void zgeev(char JOBVL,  char JOBVR, int N, std::complex<double>* A,  int LDA, std::complex<double>* W, std::complex<double>* VL, int LDVL, std::complex<double>* VR, int LDVR, std::complex<double>* WORK, int LWORK, double* RWORK);


     /**************************
      ***   HEEVD-routines
      **************************/

     void ssyevd( char JOBZ,  char UPLO,  int N, float*  Matrix,  int LDA, float* eigenvalues , float* WORK , int LWORK, int* IWORK, int LIWORK);
     void dsyevd( char JOBZ,  char UPLO,  int N, double* Matrix,  int LDA, double* eigenvalues, double* WORK, int LWORK, int* IWORK, int LIWORK);

     void cheevd( char jobz,  char uplo,  int n, std::complex<float>*  a,  int lda, float*  w, std::complex<float>*  work,  int lwork, float*  rwork, int LRWORK, int* IWORK, int LIWORK);
     void zheevd( char jobz,  char uplo,  int n, std::complex<double>* a,  int lda, double* w, std::complex<double>* work,  int lwork, double* rwork, int LRWORK, int* IWORK, int LIWORK);

     /**************************
      ***   HEEVR-routines
      **************************/

     int ssyevdx(char JOBZ,  char RANGE, char UPLO,  int N, float * A, int LDA, 
		 float VL, float VU, int IL, int UL,
		 float* W, float* WORK,  int LWORK, int* IWORK, int LIWORK);
     
     int dsyevdx(char JOBZ,  char RANGE, char UPLO,  int N, double* A, int LDA, 
		 double VL, double VU, int IL, int UL,
		 double* W, double* WORK,  int LWORK, int* IWORK, int LIWORK);
     
     int cheevdx(char JOBZ,  char RANGE, char UPLO,  int N, std::complex<float>* A, int LDA, 
		 float VL, float VU, int IL, int UL,
		 float* W, std::complex<float>* WORK, int LWORK, float* RWORK,  int LRWORK, int* IWORK, int LIWORK);
     
     int zheevdx(char JOBZ,  char RANGE, char UPLO,  int N, std::complex<double>* A, int LDA, 
		 double VL, double VU, int IL, int UL,
		 double* W, std::complex<double>* WORK,  int LWORK, double* RWORK, int LRWORK, int* IWORK, int LIWORK);

   }

   template<>
   class GEEV<GPU>
   {
   public:

     /**************************
      ***   general matrix
      **************************/

     template<typename scalartype>
     static void execute(char JOBVL, char JOBVR, 
			 matrix<scalartype, CPU>& A, 
			 vector<scalartype, CPU>& lambda_re,
			 vector<scalartype, CPU>& lambda_im,
			 matrix<scalartype, CPU>& VL,
			 matrix<scalartype, CPU>& VR);

     template<typename scalartype>
     static void execute(char JOBVL, char JOBVR, 
			 matrix<std::complex<scalartype>, CPU>& A, 
			 vector<std::complex<scalartype>, CPU>& lambda,
			 matrix<std::complex<scalartype>, CPU>& VL,
			 matrix<std::complex<scalartype>, CPU>& VR);

     /**************************
      ***   hermitian matrix
      **************************/

     template<typename scalartype>
     static void execute(char JOBZ, char UPLO, 
			 matrix<scalartype, CPU>& A, 
			 vector<scalartype, CPU>& lambda_re,
			 matrix<scalartype, CPU>& VR);

     template<typename scalartype>
     static void execute(char JOBZ, char UPLO, 
			 matrix<std::complex<scalartype>, CPU>& A, 
			 vector<             scalartype , CPU>& lambda,
			 matrix<std::complex<scalartype>, CPU>& VR);

     /*************************************
      ***   hermitian matrix (specialized)
      *************************************/

     template<typename scalartype>
     static int execute(char JOBVL, char JOBVR, 
			 scalartype LB, 
			 scalartype UB, 
			 matrix<scalartype, CPU>& A, 
			 vector<scalartype, CPU>& lambda,
			 matrix<scalartype, CPU>& VR);

     template<typename scalartype>
     static int execute(char JOBVL, char JOBVR, 
			 scalartype LB, 
			 scalartype UB, 
			 matrix<std::complex<scalartype>, CPU>& A, 
			 vector<             scalartype , CPU>& lambda,
			 matrix<std::complex<scalartype>, CPU>& V);

     template<typename scalartype>
     static int execute(char JOBVL, char JOBVR, 
			 int LB, 
			 int UB, 
			 matrix<scalartype, CPU>& A, 
			 vector<scalartype, CPU>& lambda,
			 matrix<scalartype, CPU>& VR);

     template<typename scalartype>
     static int execute(char JOBVL, char JOBVR, 
			int LB, 
			int UB, 
			matrix<std::complex<scalartype>, CPU>& A, 
			vector<             scalartype , CPU>& lambda,
			matrix<std::complex<scalartype>, CPU>& V);

   private:

     template<typename scalartype>
     static int execute(char JOBVL, char RANGE, char UPLO,
			 scalartype VL, 
			 scalartype VU, 
			 int IL, 
			 int IU, 
			 matrix<scalartype, CPU>& A, 
			 vector<scalartype, CPU>& lambda,
			 matrix<scalartype, CPU>& V);

     template<typename scalartype>
     static int execute(char JOBZ, char RANGE, char UPLO,
			 scalartype VL, 
			 scalartype VU, 
			 int IL, 
			 int IU, 
			 matrix<std::complex<scalartype>, CPU>& A, 
			 vector<             scalartype , CPU>& lambda,
			 matrix<std::complex<scalartype>, CPU>& V);

   private:

     /**************************
      ***
      ***   GEEV GPU_KERNELS_GEEV-routines
      ***
      **************************/

     static void execute(char JOBVL, char JOBVR, int N, float* A,  int LDA, float* WR, float* WI, float* VL, int LDVL, float* VR, int LDVR, float* WORK, int LWORK)
     {
       GPU_KERNELS_GEEV::sgeev(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK);
     }

     static void execute(char JOBVL, char JOBVR, int N, double* A,  int LDA, double* WR, double* WI, double* VL, int LDVL, double* VR, int LDVR, double* WORK, int LWORK )
     {
       GPU_KERNELS_GEEV::dgeev(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK);
     }

     static void execute(char JOBVL, char JOBVR, int N, std::complex<float>* A, int LDA, std::complex<float>* W, 
		  std::complex<float>* VL, int LDVL, std::complex<float>* VR, int LDVR, 
		  std::complex<float>* WORK, int LWORK, float* RWORK )
     {
       GPU_KERNELS_GEEV::cgeev(JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK);
     }

     static void execute(char JOBVL, char JOBVR, int N, std::complex<double>* A, int LDA, std::complex<double>* W, 
			   std::complex<double>* VL, int LDVL, std::complex<double>* VR, int LDVR, 
			   std::complex<double>* WORK, int LWORK, double* RWORK)
     {
       GPU_KERNELS_GEEV::zgeev(JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK);
     }

   /**************************
    ***
    ***   HEEVD GPU_KERNELS_GEEV-routines
    ***
    **************************/

     static void execute(char JOBZ,  char UPLO,  int N, float*  A,  int LDA, float*  W, 
			 float*  WORK, int LWORK, int* IWORK, int LIWORK)
     {
       GPU_KERNELS_GEEV::ssyevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK);
     }

     static void execute(char JOBZ,  char UPLO,  int N, double* A,  int LDA, double* W, 
		  double* WORK, int LWORK, int* IWORK, int LIWORK)
     {
       GPU_KERNELS_GEEV::dsyevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK);
     }

     static void execute(char JOBZ,  char UPLO, int N, std::complex<float>*  A,  int LDA, float*  W, 
		  std::complex<float>* WORK , int LWORK, 
		  float*               RWORK, int LRWORK,
		  int*                 IWORK, int LIWORK)
     {
       GPU_KERNELS_GEEV::cheevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK);
     }

     static void execute(char JOBZ,  char UPLO, int N, std::complex<double>* A,  int LDA, double* W, 
		  std::complex<double>* WORK , int LWORK, 
		  double*               RWORK, int LRWORK,
		  int*                  IWORK, int LIWORK)
     {
       GPU_KERNELS_GEEV::zheevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK);
     }

   /****************************
    ***
    ***   HEEVDX GPU_KERNELS_GEEV-routines
    ***
    ****************************/

     static int execute(char JOBZ,  char RANGE, char UPLO,  int N, float* A, int LDA, 
			float VL, float VU, int IL, int UL,
			float* W, float* WORK,  int LWORK, int* IWORK, int LIWORK)
     {
       return GPU_KERNELS_GEEV::ssyevdx(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, UL,
				       W, WORK, LWORK, IWORK, LIWORK);
     }
     
     static int execute(char JOBZ,  char RANGE, char UPLO,  int N, double* A, int LDA, 
			double VL, double VU, int IL, int UL,
			double* W, double* WORK, int LWORK, int* IWORK, int LIWORK)
     {
       return GPU_KERNELS_GEEV::dsyevdx(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, UL,
				       W, WORK, LWORK, IWORK, LIWORK);
     }
     
     static int execute(char JOBZ,  char RANGE, char UPLO,  int N, std::complex<float>* A, int LDA, 
			float VL, float VU, int IL, int UL,
			float* W, std::complex<float>* WORK, int LWORK, float* RWORK, int LRWORK, int* IWORK, int LIWORK)
     {
       return GPU_KERNELS_GEEV::cheevdx(JOBZ, RANGE, UPLO, N, A, LDA, 
				       VL, VU, IL, UL,
				       W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK);
     }
     
     static int execute(char JOBZ,  char RANGE, char UPLO,  int N, std::complex<double>* A, int LDA, 
			double VL, double VU, int IL, int UL,
			double* W, std::complex<double>* WORK, int LWORK, double* RWORK, int LRWORK, int* IWORK, int LIWORK)
     {
       return GPU_KERNELS_GEEV::zheevdx(JOBZ, RANGE, UPLO, N, A, LDA, 
				       VL, VU, IL, UL,
				       W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK);
     }
   };

   /**************************
    ***
    ***   GEEV-routines
    ***
    **************************/

   template<typename scalartype>
   void GEEV<GPU>::execute(char JOBVL, char JOBVR, 
			   matrix<scalartype, CPU>& A, 
			   vector<scalartype, CPU>& lambda_re,
			   vector<scalartype, CPU>& lambda_im,
			   matrix<scalartype, CPU>& VL,
			   matrix<scalartype, CPU>& VR)
   {
     if( JOBVL != 'N' and JOBVL != 'V')
       throw std::logic_error(__FUNCTION__);

     if( JOBVR != 'N' and JOBVR != 'V')
       throw std::logic_error(__FUNCTION__);

     if( A.get_current_size().first !=  A.get_current_size().second)
       throw std::logic_error(__FUNCTION__);

     matrix<std::complex<scalartype>, CPU> X;
     X.copy_from(A);

     int N_A = A.get_current_size().first;
     int LDA = A.get_global_size().first;

     int LDVL = VL.get_global_size().first;
     int LDVR = VR.get_global_size().first;

     int LWORK = -1;
     vector<scalartype, CPU> WORK(1);

     {
       execute(JOBVL, JOBVR, N_A, &X(0,0), LDA, &lambda_re(0), &lambda_im(0), &VL(0,0), LDVL, &VR(0,0), LDVR, &WORK, LWORK);

       LWORK = WORK[0];
       WORK.resize(LWORK);
     }

     execute(JOBVL, JOBVR, N_A, &X(0,0), LDA, &lambda_re(0), &lambda_im(0), &VL(0,0), LDVL, &VR(0,0), LDVR, &WORK[0], LWORK);
   }

   template<typename scalartype>
   void GEEV<GPU>::execute(char JOBVL, char JOBVR, 
			   matrix<std::complex<scalartype>, CPU>& A, 
			   vector<std::complex<scalartype>, CPU>& lambda,
			   matrix<std::complex<scalartype>, CPU>& VL,
			   matrix<std::complex<scalartype>, CPU>& VR)
   {
     if( JOBVL != 'N' and JOBVL != 'V')
       throw std::logic_error(__FUNCTION__);

     if( JOBVR != 'N' and JOBVR != 'V')
       throw std::logic_error(__FUNCTION__);

     if( A.get_current_size().first !=  A.get_current_size().second)
       throw std::logic_error(__FUNCTION__);

     matrix<std::complex<scalartype>, CPU> X;
     X.copy_from(A);

     int N_A = A.get_current_size().first;
     int LDA = A.get_global_size().first;

     int LDVL = VL.get_global_size().first;
     int LDVR = VR.get_global_size().first;

     int LWORK = -1;
     vector<std::complex<scalartype>, CPU> WORK(1);

     vector<std::complex<scalartype>, CPU> RWORK(2*N_A);

     {
       std::complex<scalartype> WORK;
       execute(JOBVL, JOBVR, N_A, &X(0,0), LDA, &lambda(0), &VL(0,0), LDVL, &VR(0,0), LDVR, &WORK, LWORK, &RWORK(0));

       LWORK = real(WORK[0]);
       WORK.resize(LWORK);
     }

     execute(JOBVL, JOBVR, N_A, &X(0,0), LDA, &lambda(0), &VL(0,0), LDVL, &VR(0,0), LDVR, &WORK, LWORK, &RWORK(0));
   }

   /**************************
    ***
    ***   HEEVD-routines
    ***
    **************************/

     template<typename scalartype>
     void GEEV<GPU>::execute(char JOBZ, char UPLO, 
			     matrix<scalartype, CPU>& A, 
			     vector<scalartype, CPU>& lambda,
			     matrix<scalartype, CPU>& VR)
     {
       if( JOBZ != 'N' and JOBZ != 'V')
	 throw std::logic_error(__FUNCTION__);

       if( UPLO != 'U' and UPLO != 'L')
	 throw std::logic_error(__FUNCTION__);

       int N_A = A.get_current_size().first;
       int LDA = A.get_global_size().first;

       VR.copy_from(A);

       int LWORK  = -1;
       vector<scalartype, CPU> WORK (1);

       int LIWORK = -1;
       vector<int , CPU> IWORK("IWORK", 1);

       {
	 execute(JOBZ, UPLO, N_A, &VR(0,0), LDA, &lambda(0), &WORK, LWORK, &IWORK[0], LIWORK);

	 LWORK = real(WORK[0]);
	 WORK.resize(LWORK);

	 LIWORK = IWORK[0];
	 IWORK.resize(LIWORK);
       }
       
       execute(JOBZ, UPLO, N_A, &VR(0,0), LDA, &lambda(0), &WORK[0], LWORK, &IWORK[0], LIWORK);
     }

     template<typename scalartype>
     void GEEV<GPU>::execute(char JOBZ, char UPLO, 
			     matrix<std::complex<scalartype>, CPU>& A, 
			     vector<             scalartype , CPU>& lambda,
			     matrix<std::complex<scalartype>, CPU>& VR)
     {
       if( JOBZ != 'N' and JOBZ != 'V')
	 throw std::logic_error(__FUNCTION__);

       if( UPLO != 'U' and UPLO != 'L')
	 throw std::logic_error(__FUNCTION__);

       VR.copy_from(A);

       int N_A = VR.get_current_size().first;
       int LDA = VR.get_global_size().first;

       int LWORK = -1;
       vector<std::complex<scalartype>, CPU> WORK("WORK", 1);

       int LRWORK = -1;
       vector<scalartype , CPU> RWORK("RWORK", 1);

       int LIWORK = -1;
       vector<int , CPU> IWORK("IWORK", 1);
       
       {
	 execute(JOBZ, UPLO, N_A, &VR(0,0), LDA, &lambda[0], &WORK[0], LWORK, &RWORK[0], LRWORK, &IWORK[0], LIWORK);
	 
	 LWORK = real(WORK[0]);
	 WORK.resize(LWORK);
	 
	 LRWORK = RWORK[0];
	 RWORK.resize(LRWORK);
	 
	 LIWORK = IWORK[0];
	 IWORK.resize(LIWORK);
       }
       
       execute(JOBZ, UPLO, N_A, &VR(0,0), LDA, &lambda[0], &WORK[0], LWORK, &RWORK[0], LRWORK, &IWORK[0], LIWORK);
     }

   /**************************
    ***
    ***   HEEVDX-routines
    ***
    **************************/

   template<typename scalartype>
   int GEEV<GPU>::execute(char JOBZ, char UPLO,
			  scalartype VL, 
			  scalartype VU, 
			  matrix<scalartype, CPU>& A, 
			  vector<scalartype, CPU>& lambda,
			  matrix<scalartype, CPU>& VR)
   {
     char RANGE = 'V';
     
     int IL = -1;
     int IU = -1;

     return execute(JOBZ, RANGE, UPLO, VL, VU, IL, IU, A, lambda, VR);
   }

   template<typename scalartype>
   int GEEV<GPU>::execute(char JOBZ, char UPLO,
			  scalartype VL, 
			  scalartype VU, 
			  matrix<std::complex<scalartype>, CPU>& A, 
			  vector<             scalartype , CPU>& lambda,
			  matrix<std::complex<scalartype>, CPU>& VR)
   {
     char RANGE = 'V';

     int IL = -1;
     int IU = -1;

     return execute(JOBZ, RANGE, UPLO, VL, VU, IL, IU, A, lambda, VR);
   }

   template<typename scalartype>
   int GEEV<GPU>::execute(char JOBZ, char UPLO,
			  int IL, int IU, 
			  matrix<scalartype, CPU>& A, 
			  vector<scalartype, CPU>& lambda,
			  matrix<scalartype, CPU>& VR)
   {    
     char RANGE = 'I';
     
     scalartype VL = 0.;
     scalartype VU = 0.;

     return execute(JOBZ, RANGE, UPLO, VL, VU, IL, IU,  A, lambda, VR);
   }

   template<typename scalartype>
   int GEEV<GPU>::execute(char JOBZ, char UPLO,
			  int IL, int IU, 
			  matrix<std::complex<scalartype>, CPU>& A, 
			  vector<             scalartype , CPU>& lambda,
			  matrix<std::complex<scalartype>, CPU>& VR)
   {    
     char RANGE = 'I';
     
     scalartype VL = 0.;
     scalartype VU = 0.;

     return execute(JOBZ, RANGE, UPLO, VL, VU, IL, IU,  A, lambda, VR);
   }

   template<typename scalartype>
   int GEEV<GPU>::execute(char JOBZ, char RANGE, char UPLO,
			  scalartype VL, 
			  scalartype VU, 
			  int IL, 
			  int IU, 
			  matrix<scalartype, CPU>& A, 
			  vector<scalartype, CPU>& lambda,
			  matrix<scalartype, CPU>& V)
   {
     if( JOBZ != 'N' and JOBZ != 'V')
       throw std::logic_error(__FUNCTION__);

     if( UPLO != 'U' and UPLO != 'L')
       throw std::logic_error(__FUNCTION__);

     V.copy_from(A);

     int N   = V.get_current_size().first;
     int LDA = V.get_global_size().first;

     int M=-1;

     int LWORK  = -1;
     int LIWORK = -1;

     vector<scalartype, CPU> WORK (1);
     vector<int       , CPU> IWORK(1);

     {// find optimal work-space
       execute(JOBZ, RANGE, UPLO, N, &V(0,0), LDA, VL, VU, IL, IU, &lambda[0], &WORK[0], LWORK, &IWORK[0], LIWORK);
       
       LWORK = WORK[0];       
       WORK.resize(LWORK);

       LIWORK = IWORK[0];       
       IWORK.resize(LIWORK);
     }
     
     M = execute(JOBZ, RANGE, UPLO, N, &V(0,0), LDA, VL, VU, IL, IU, &lambda[0], &WORK[0], LWORK, &IWORK[0], LIWORK);

     return M;
   }

   template<typename scalartype>
   int GEEV<GPU>::execute(char JOBZ, char RANGE, char UPLO,
			   scalartype VL, 
			   scalartype VU, 
			   int IL, 
			   int IU, 
			   matrix<std::complex<scalartype>, CPU>& A, 
			   vector<             scalartype , CPU>& lambda,
			   matrix<std::complex<scalartype>, CPU>& V)
   {
     if( JOBZ != 'N' and JOBZ != 'V')
       throw std::logic_error(__FUNCTION__);

     if( UPLO != 'U' and UPLO != 'L')
       throw std::logic_error(__FUNCTION__);

     V.copy_from(A);

     int N   = V.get_current_size().first;
     int LDA = V.get_global_size().first;

     int M      = -1;

     int LWORK  = -1;
     int LRWORK = -1;
     int LIWORK = -1;

     vector<std::complex<scalartype>, CPU> WORK (1);
     vector<             scalartype , CPU> RWORK(1);
     vector<             int        , CPU> IWORK(1);
     
     {// find optimal work-space
       execute(JOBZ, RANGE, UPLO, N, &V(0,0), LDA, VL, VU, IL, IU, &lambda[0], &WORK[0], LWORK, &RWORK[0], LRWORK, &IWORK[0], LIWORK);
       
       LWORK = real(WORK[0]);       
       WORK.resize(LWORK);

       LRWORK = RWORK[0];       
       RWORK.resize(LRWORK);

       LIWORK = IWORK[0];       
       IWORK.resize(LIWORK);
     }
     
     M = execute(JOBZ, RANGE, UPLO, N, &V(0,0), LDA, VL, VU, IL, IU, &lambda[0], &WORK[0], LWORK, &RWORK[0], LRWORK, &IWORK[0], LIWORK);

     return M;
   }

 }

#endif
