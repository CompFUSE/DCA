//-*-C++-*-

#ifndef LINALG_GEEV_CPU_H
#define LINALG_GEEV_CPU_H

namespace LIN_ALG {

  template<>
  class GEEV<CPU>
  {
  public:

    /**************************
     ***   general matrix
     **************************/

    template<typename scalartype>
    static void execute(char JOBVL, char JOBVR,
                        matrix<scalartype, CPU>& A,
                        dca::linalg::Vector<scalartype, CPU>& lambda_re,
                        dca::linalg::Vector<scalartype, CPU>& lambda_im,
                        matrix<scalartype, CPU>& VL,
                        matrix<scalartype, CPU>& VR);

    template<typename scalartype>
    static void execute(char JOBVL, char JOBVR,
                        matrix<std::complex<scalartype>, CPU>& A,
                        dca::linalg::Vector<std::complex<scalartype>, CPU>& lambda,
                        matrix<std::complex<scalartype>, CPU>& VL,
                        matrix<std::complex<scalartype>, CPU>& VR);

    /**************************
     ***   hermitian matrix
     **************************/

    template<typename scalartype>
    static void execute(char JOBZ, char UPLO,
                        matrix<scalartype, CPU>& A,
                        dca::linalg::Vector<scalartype, CPU>& lambda_re,
                        matrix<scalartype, CPU>& VR);

    template<typename scalartype>
    static void execute(char JOBZ, char UPLO,
                        matrix<std::complex<scalartype>, CPU>& A,
                        dca::linalg::Vector<scalartype, CPU>& lambda,
                        matrix<std::complex<scalartype>, CPU>& VR);

    template<typename scalartype>
    static void execute_on_Greens_function_matrix(char JOBZ, char UPLO,
                                                  matrix<std::complex<scalartype>, CPU>& A,
                                                  dca::linalg::Vector<scalartype, CPU>& lambda,
                                                  matrix<std::complex<scalartype>, CPU>& VR);

    /*
    template<typename scalartype>
    class hermitian_plan
    {
    public:

      hermitian_plan();
      ~hermitian_plan();

      void initialize(char JOBZ, char UPLO,
		      matrix<scalartype, CPU>& A,
		      dca::linalg::Vector<scalartype, CPU>& lambda_re,
		      matrix<scalartype, CPU>& VR);

      void execute(char JOBZ, char UPLO,
		   matrix<scalartype, CPU>& A,
		   dca::linalg::Vector<scalartype, CPU>& lambda_re,
		   matrix<scalartype, CPU>& VR);

    private:

      bool initialzed;

      int LIWORK;
      dca::linalg::Vector<int, CPU> IWORK;

      int LWORK;
      dca::linalg::Vector<scalartype, CPU> WORK;      
    };
    */

    template<typename scalartype>
    class hermitian_plan
    {
    public:

      hermitian_plan();
      ~hermitian_plan();

      void initialize(char JOBZ, char UPLO,
		      matrix<scalartype, CPU>& A,
		      dca::linalg::Vector<scalartype, CPU>& lambda_re,
		      matrix<scalartype, CPU>& VR);

      void initialize(char JOBZ, char UPLO,
		      matrix<std::complex<scalartype>, CPU>& A,
		      dca::linalg::Vector<std::complex<scalartype>, CPU>& lambda_re,
		      matrix<std::complex<scalartype>, CPU>& VR);

      void execute(char JOBZ, char UPLO,
		   matrix<scalartype, CPU>& A,
		   dca::linalg::Vector<scalartype, CPU>& lambda_re,
		   matrix<scalartype, CPU>& VR);

      void execute(char JOBZ, char UPLO,
		   matrix<std::complex<scalartype>, CPU>& A,
		   dca::linalg::Vector<std::complex<scalartype>, CPU>& lambda_re,
		   matrix<std::complex<scalartype>, CPU>& VR);

    private:

      bool initialzed;

      int LIWORK;
      dca::linalg::Vector<int, CPU> IWORK;

      int LRWORK;
      dca::linalg::Vector<scalartype, CPU> RWORK;      

      int LWORK;
      dca::linalg::Vector<std::complex<scalartype>, CPU> WORK;      
    };

    /*************************************
     ***   hermitian matrix (specialized)
     *************************************/

    template<typename scalartype>
    static int execute(char JOBVL, char JOBVR,
                       scalartype LB,
                       scalartype UB,
                       matrix<scalartype, CPU>& A,
                       dca::linalg::Vector<scalartype, CPU>& lambda,
                       matrix<scalartype, CPU>& VR,
                       std::string              alg_type="DEFAULT");

    template<typename scalartype>
    static int execute(char JOBVL, char JOBVR,
                       scalartype LB,
                       scalartype UB,
                       matrix<std::complex<scalartype>, CPU>& A,
                       dca::linalg::Vector<scalartype, CPU>& lambda,
                       matrix<std::complex<scalartype>, CPU>& V,
                       std::string              alg_type="DEFAULT");

    template<typename scalartype>
    static int execute(char JOBVL, char JOBVR,
                       int LB,
                       int UB,
                       matrix<scalartype, CPU>& A,
                       dca::linalg::Vector<scalartype, CPU>& lambda,
                       matrix<scalartype, CPU>& VR,
                       std::string              alg_type="DEFAULT");

    template<typename scalartype>
    static int execute(char JOBVL, char JOBVR,
                       int LB,
                       int UB,
                       matrix<std::complex<scalartype>, CPU>& A,
                       dca::linalg::Vector<scalartype, CPU>& lambda,
                       matrix<std::complex<scalartype>, CPU>& V,
                       std::string              alg_type="DEFAULT");

    /**************************
     ***   hermitian matrix
     **************************/

    template<typename scalartype>
    static void execute_on_small_matrix(char JOBZ, char UPLO,
                                        matrix<scalartype, CPU>& A,
                                        dca::linalg::Vector<scalartype, CPU>& lambda_re,
                                        matrix<scalartype, CPU>& VR);

    template<typename scalartype>
    static void execute_on_small_matrix(char JOBZ, char UPLO,
                                        matrix<std::complex<scalartype>, CPU>& A,
                                        dca::linalg::Vector<scalartype, CPU>& lambda,
                                        matrix<std::complex<scalartype>, CPU>& VR);

  private:

    template<typename scalartype>
    static int execute(char JOBVL, char RANGE, char UPLO,
                       scalartype VL,
                       scalartype VU,
                       int IL,
                       int IU,
                       matrix<scalartype, CPU>& A,
                       dca::linalg::Vector<scalartype, CPU>& lambda,
                       matrix<scalartype, CPU>& V,
                       std::string              alg_type);

    template<typename scalartype>
    static int execute(char JOBZ, char RANGE, char UPLO,
                       scalartype VL,
                       scalartype VU,
                       int IL,
                       int IU,
                       matrix<std::complex<scalartype>, CPU>& A,
                       dca::linalg::Vector<scalartype, CPU>& lambda,
                       matrix<std::complex<scalartype>, CPU>& V,
                       std::string              alg_type);

  public:

    /**************************
     ***   GEEV-routines
     **************************/

    static void execute(char JOBVL, char JOBVR, int N, float*  A, int LDA, float*   WR, float*   WI, float* VL , int LDVL, float* VR , int LDVR, float* WORK , int LWORK, int INFO);
    static void execute(char JOBVL, char JOBVR, int N, double* A, int LDA, double*  WR, double*  WI, double* VL, int LDVL, double* VR, int LDVR, double* WORK, int LWORK, int INFO);

    static void execute(char JOBVL, char JOBVR, int N, std::complex<float>*  A, int LDA, std::complex<float>* W,
                        std::complex<float>*  VL, int LDVL, std::complex<float>*  VR, int LDVR,
                        std::complex<float>*  WORK, int LWORK, float*  RWORK, int INFO );

    static void execute(char JOBVL, char JOBVR, int N, std::complex<double>* A, int LDA, std::complex<double>* W,
                        std::complex<double>* VL, int LDVL, std::complex<double>* VR, int LDVR,
                        std::complex<double>* WORK, int LWORK, double* RWORK, int INFO );

    /**************************
     ***   HEEV-routines
     **************************/

    static void execute(char JOBZ,  char UPLO,  int N, float*  A,  int LDA, float*  W, float*  WORK,  int LWORK);
    static void execute(char JOBZ,  char UPLO,  int N, double* A,  int LDA, double* W, double* WORK,  int LWORK);

    static void execute(char jobz,  char uplo,  int n, std::complex<float>*  A,  int lda, float*  w, std::complex<float>*  work,  int lwork, float*  rwork );
    static void execute(char jobz,  char uplo,  int n, std::complex<double>* A,  int lda, double* w, std::complex<double>* work,  int lwork, double* rwork );

    /**************************
     ***   HEEVD-routines
     **************************/

    static void execute(char JOBZ,  char UPLO,  int N, float*  A,  int LDA, float*  W,
                        float*  WORK, int LWORK, int* IWORK, int LIWORK);

    static void execute(char JOBZ,  char UPLO,  int N, double* A,  int LDA, double* W,
                        double* WORK, int LWORK, int* IWORK, int LIWORK);

    static void execute(char jobz,  char uplo, int n, std::complex<float>*  A,  int lda, float*  w,
                        std::complex<float>* work , int lwork,
                        float*               rwork, int lrwork,
                        int*                 iwork, int liwork);

    static void execute(char jobz,  char uplo, int n, std::complex<double>* A,  int lda, double* w,
                        std::complex<double>* work , int lwork,
                        double*               rwork, int lrwork,
                        int*                  iwork, int liwork);

    /**************************
     ***   HEEVX-routines
     **************************/

    static int execute(char JOBZ,  char RANGE, char UPLO,  int N, float * A, int LDA,
                       float VL, float VU, int IL, int UL,
                       float* W, float* Z, int LDZ,
                       float* WORK, int LWORK, int* IWORK, int* IFAIL, int INFO );

    static int execute(char JOBZ,  char RANGE, char UPLO,  int N, double* A, int LDA,
                       double VL, double VU, int IL, int UL,
                       double* W, double* Z, int LDZ,
                       double* WORK, int LWORK, int* IWORK, int* IFAIL, int INFO );

    static int execute(char JOBZ,  char RANGE, char UPLO, int N, std::complex<float>* A, int LDA,
                       float VL, float VU, int IL, int UL,
                       float* W, std::complex<float>* Z, int LDZ,
                       std::complex<float>* WORK, int LWORK, float* RWORK, int LRWORK, int* IWORK, int* IFAIL, int INFO );

    static int execute(char JOBZ,  char RANGE, char UPLO, int N, std::complex<double>* A, int LDA,
                       double VL, double VU, int IL, int UL,
                       double* W, std::complex<double>* Z, int LDZ,
                       std::complex<double>* WORK, int LWORK, double* RWORK, int LRWORK, int* IWORK, int* IFAIL, int INFO );

    /**************************
     ***   HEEVR-routines
     **************************/

    static void execute(char JOBZ,  char RANGE, char UPLO,  int N, float* A, int LDA,
                        float VL, float VU, int IL, int UL, int M,
                        float* W, float* Z, int LDZ, int* ISUPPZ,
                        float* WORK,  int LWORK, int* IWORK, int LIWORK, int INFO );

    static void execute(char JOBZ,  char RANGE, char UPLO,  int N, double* A, int LDA,
                        double VL, double VU, int IL, int UL, int M,
                        double* W, double* Z, int LDZ, int* ISUPPZ,
                        double* WORK, int LWORK, int* IWORK, int LIWORK, int INFO );

    static void execute(char JOBZ,  char RANGE, char UPLO,  int N, std::complex<float>* A, int LDA,
                        float VL, float VU, int IL, int UL, int M,
                        float* W, std::complex<float>* Z, int LDZ, int* ISUPPZ,
                        std::complex<float>* WORK, int LWORK, float* RWORK, int LRWORK, int* IWORK, int LIWORK, int INFO );

    static void execute(char JOBZ,  char RANGE, char UPLO,  int N, std::complex<double>* A, int LDA,
                        double VL, double VU, int IL, int UL, int M,
                        double* W, std::complex<double>* Z, int LDZ, int* ISUPPZ,
                        std::complex<double>* WORK, int LWORK, double* RWORK, int LRWORK, int* IWORK, int LIWORK, int INFO );
  };

  /**************************
   ***
   ***   GEEV-routines
   ***
   **************************/

  template<typename scalartype>
  void GEEV<CPU>::execute(char JOBVL, char JOBVR,
                          matrix<scalartype, CPU>& A,
                          dca::linalg::Vector<scalartype, CPU>& lambda_re,
                          dca::linalg::Vector<scalartype, CPU>& lambda_im,
                          matrix<scalartype, CPU>& VL,
                          matrix<scalartype, CPU>& VR)
  {
    if( JOBVL != 'N' and JOBVL != 'V')
      throw std::logic_error(__FUNCTION__);

    if( JOBVR != 'N' and JOBVR != 'V')
      throw std::logic_error(__FUNCTION__);

    if( A.size().first !=  A.size().second)
      throw std::logic_error(__FUNCTION__);

    matrix<std::complex<scalartype>, CPU> X(A);

    int N_A = A.size().first;
    int LDA = A.leadingDimension();

    int LDVL = VL.leadingDimension();
    int LDVR = VR.leadingDimension();

    int LWORK = -1;

    int INFO = -1;

    {
      scalartype WORK;
      execute(JOBVL, JOBVR, N_A, &X(0,0), LDA, &lambda_re(0), &lambda_im(0), &VL(0,0), LDVL, &VR(0,0), LDVR, &WORK, LWORK, INFO);

      LWORK = WORK[0];
    }

    dca::linalg::Vector<scalartype, CPU> WORK(LWORK);

    execute(JOBVL, JOBVR, N_A, &X(0,0), LDA, &lambda_re(0), &lambda_im(0), &VL(0,0), LDVL, &VR(0,0), LDVR, &WORK[0], LWORK, INFO);
  }

  template<typename scalartype>
  void GEEV<CPU>::execute(char JOBVL, char JOBVR,
                          matrix<std::complex<scalartype>, CPU>& A,
                          dca::linalg::Vector<std::complex<scalartype>, CPU>& lambda,
                          matrix<std::complex<scalartype>, CPU>& VL,
                          matrix<std::complex<scalartype>, CPU>& VR)
  {
    if( JOBVL != 'N' and JOBVL != 'V')
      throw std::logic_error(__FUNCTION__);

    if( JOBVR != 'N' and JOBVR != 'V')
      throw std::logic_error(__FUNCTION__);

    if( A.size().first !=  A.size().second)
      throw std::logic_error(__FUNCTION__);

    matrix<std::complex<scalartype>, CPU> X(A);

    int N_A = A.size().first;
    int LDA = A.leadingDimension();

    int LDVL = VL.leadingDimension();
    int LDVR = VR.leadingDimension();

    int LWORK = -1;

    dca::linalg::Vector<scalartype, CPU> RWORK(2*N_A);

    int INFO = -1;

    {
      std::complex<scalartype> WORK;
      execute(JOBVL, JOBVR, N_A, &X(0,0), LDA, &lambda[0], &VL(0,0), LDVL, &VR(0,0), LDVR, &WORK, LWORK, &RWORK[0], INFO);

      LWORK = real(WORK);
    }

    dca::linalg::Vector<std::complex<scalartype>, CPU> WORK(LWORK);

    execute(JOBVL, JOBVR, N_A, &X(0,0), LDA, &lambda[0], &VL(0,0), LDVL, &VR(0,0), LDVR, &WORK[0], LWORK, &RWORK[0], INFO);
  }

  /**************************
   ***
   ***   HEEV-routines
   ***
   **************************/

  template<typename scalartype>
  void GEEV<CPU>::execute(char JOBZ, char UPLO,
                          matrix<scalartype, CPU>& A,
                          dca::linalg::Vector<scalartype, CPU>& lambda,
                          matrix<scalartype, CPU>& VR)
  {
    if( JOBZ != 'N' and JOBZ != 'V')
      throw std::logic_error(__FUNCTION__);

    if( UPLO != 'U' and UPLO != 'L')
      throw std::logic_error(__FUNCTION__);

    int N_A = A.size().first;
    int LDA = A.leadingDimension();

    VR = A;

    if(false)
      {
        int LWORK = -1;
        dca::linalg::Vector<scalartype, CPU> WORK("WORK", 1);

        {
          execute(JOBZ, UPLO, N_A, &VR(0,0), LDA, &lambda[0], &WORK[0], LWORK);

          LWORK = WORK[0]*(1.+1.e-3);
          WORK.resize(LWORK);
        }

        execute(JOBZ, UPLO, N_A, &VR(0,0), LDA, &lambda[0], &WORK[0], LWORK);
      }
    else
      {
        int LWORK = -1;
        dca::linalg::Vector<scalartype, CPU> WORK("WORK", 1);

        int LIWORK = -1;
        dca::linalg::Vector<int, CPU> IWORK("LIWORK", 1);

        {
          execute(JOBZ, UPLO, N_A, &VR(0,0), LDA, &lambda[0], &WORK[0], LWORK, &IWORK[0], LIWORK);

          LWORK = WORK[0]*(1.+1.e-3);
          WORK.resize(LWORK);

          LIWORK = IWORK[0]*(1.+1.e-3);
          IWORK.resize(LIWORK);
        }

        execute(JOBZ, UPLO, N_A, &VR(0,0), LDA, &lambda[0], &WORK[0], LWORK, &IWORK[0], LIWORK);
      }
  }

  template<typename scalartype>
  void GEEV<CPU>::execute(char JOBZ, char UPLO,
                          matrix<std::complex<scalartype>, CPU>& A,
                          dca::linalg::Vector<scalartype, CPU>& lambda,
                          matrix<std::complex<scalartype>, CPU>& VR)
  {
    if( JOBZ != 'N' and JOBZ != 'V')
      throw std::logic_error(__FUNCTION__);

    if( UPLO != 'U' and UPLO != 'L')
      throw std::logic_error(__FUNCTION__);

    VR = A;

    int N_A = VR.size().first;
    int LDA = VR.leadingDimension();

    if(false)
      {
        int LWORK = -1;
        dca::linalg::Vector<std::complex<scalartype>, CPU> WORK("WORK", 1);

        int LRWORK = std::max(1, 3*N_A-2);
        dca::linalg::Vector<scalartype , CPU> RWORK("RWORK", LRWORK);

        {
          execute(JOBZ, UPLO, N_A, &VR(0,0), LDA, &lambda[0], &WORK[0], LWORK, &RWORK[0]);

          LWORK = real(WORK[0]);

          WORK.resize(LWORK);
        }

        execute(JOBZ, UPLO, N_A, &VR(0,0), LDA, &lambda[0], &WORK[0], LWORK, &RWORK[0]);
      }
    else
      {
        int LWORK = -1;
        dca::linalg::Vector<std::complex<scalartype>, CPU> WORK("WORK", 1);

        int LRWORK = -1;
        dca::linalg::Vector<scalartype , CPU> RWORK("RWORK", 1);

        int LIWORK = -1;
        dca::linalg::Vector<int, CPU> IWORK("RWORK", 1);

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
  }

  template<typename scalartype>
  void GEEV<CPU>::execute_on_Greens_function_matrix(char JOBZ, char UPLO,
                                                    matrix<std::complex<scalartype>, CPU>& A,
                                                    dca::linalg::Vector<scalartype, CPU>& lambda,
                                                    matrix<std::complex<scalartype>, CPU>& VR)
  {
    int N = A.size().first;
    assert(N%2==0);

    switch(N)
      {
      case 2:
        {
	    assert(std::fabs(imag(A(0,0)))<1.e-6);
	    assert(std::fabs(imag(A(1,1)))<1.e-6);

          lambda[0] = real(A(0,0));
          lambda[1] = real(A(1,1));

          VR(0,0) = 1.0; VR(0,1) = 0.0;
          VR(1,0) = 0.0; VR(1,1) = 1.0;
        }
        break;

      default:
        execute(JOBZ, UPLO, A, lambda, VR);
      }
  }

  /**************************
   ***
   ***   hermitian matrix
   ***
   **************************/

  template<typename scalartype>
  void GEEV<CPU>::execute_on_small_matrix(char JOBZ, char UPLO,
                                          matrix<scalartype, CPU>& A,
                                          dca::linalg::Vector<scalartype, CPU>& lambda,
                                          matrix<scalartype, CPU>& VR)
  {
    int N = A.size().first;

    switch(N)
      {
      case 1:
        {
          lambda[0] = A(0,0);

          VR(0,0)   = 1.0;
        }
        break;

      default:
        execute(JOBZ, UPLO, A, lambda, VR);
      }
  }

  template<typename scalartype>
  void GEEV<CPU>::execute_on_small_matrix(char JOBZ, char UPLO,
                                          matrix<std::complex<scalartype>, CPU>& A,
                                          dca::linalg::Vector<scalartype, CPU>& lambda,
                                          matrix<std::complex<scalartype>, CPU>& VR)
  {
    int N = A.size().first;

    switch(N)
      {
      case 1:
        {
          lambda[0] = real(A(0,0));

          real(VR(0,0)) = 1.0;
          imag(VR(0,0)) = 0.0;
        }
        break;

      default:
        execute(JOBZ, UPLO, A, lambda, VR);
      }
  }


  /**************************
   ***
   ***   HEEVx-routines
   ***
   **************************/

  template<typename scalartype>
  int GEEV<CPU>::execute(char JOBZ, char UPLO,
                         scalartype VL,
                         scalartype VU,
                         matrix<scalartype, CPU>& A,
                         dca::linalg::Vector<scalartype, CPU>& lambda,
                         matrix<scalartype, CPU>& VR,
                         std::string              alg_type)
  {
    char RANGE = 'V';

    int IL = -1;
    int IU = -1;

    return execute(JOBZ, RANGE, UPLO, VL, VU, IL, IU, A, lambda, VR, alg_type);
  }

  template<typename scalartype>
  int GEEV<CPU>::execute(char JOBZ, char UPLO,
                         scalartype VL,
                         scalartype VU,
                         matrix<std::complex<scalartype>, CPU>& A,
                         dca::linalg::Vector<scalartype, CPU>& lambda,
                         matrix<std::complex<scalartype>, CPU>& VR,
                         std::string              alg_type)
  {
    char RANGE = 'V';

    int IL = -1;
    int IU = -1;

    return execute(JOBZ, RANGE, UPLO, VL, VU, IL, IU, A, lambda, VR, alg_type);
  }

  template<typename scalartype>
  int GEEV<CPU>::execute(char JOBZ, char UPLO,
                         int IL, int IU,
                         matrix<scalartype, CPU>& A,
                         dca::linalg::Vector<scalartype, CPU>& lambda,
                         matrix<scalartype, CPU>& VR,
                         std::string              alg_type)
  {
    char RANGE = 'I';

    scalartype VL = 0.;
    scalartype VU = 0.;

    return execute(JOBZ, RANGE, UPLO, VL, VU, IL, IU,  A, lambda, VR, alg_type);
  }

  template<typename scalartype>
  int GEEV<CPU>::execute(char JOBZ, char UPLO,
                         int IL, int IU,
                         matrix<std::complex<scalartype>, CPU>& A,
                         dca::linalg::Vector<scalartype, CPU>& lambda,
                         matrix<std::complex<scalartype>, CPU>& VR,
                         std::string                            alg_type)
  {
    char RANGE = 'I';

    scalartype VL = 0.;
    scalartype VU = 0.;

    return execute(JOBZ, RANGE, UPLO, VL, VU, IL, IU,  A, lambda, VR, alg_type);
  }

  template<typename scalartype>
  int GEEV<CPU>::execute(char JOBZ, char RANGE, char UPLO,
                         scalartype VL,
                         scalartype VU,
                         int IL,
                         int IU,
                         matrix<scalartype, CPU>& A,
                         dca::linalg::Vector<scalartype, CPU>& lambda,
                         matrix<scalartype, CPU>& V,
                         std::string              alg_type)
  {
    if( JOBZ != 'N' and JOBZ != 'V')
      throw std::logic_error(__FUNCTION__);

    if( UPLO != 'U' and UPLO != 'L')
      throw std::logic_error(__FUNCTION__);

    matrix<scalartype, CPU> X(A);

    V = A;

    int N_A = A.size().first;

    int LDA = X.leadingDimension();
    int LDZ = V.leadingDimension();

    int M=-1;

    if(alg_type == "DEFAULT")
      {
        dca::linalg::Vector<scalartype, CPU> WORK(1);

        dca::linalg::Vector<int, CPU> IWORK(5*N_A);
        dca::linalg::Vector<int, CPU> IFAIL(N_A);

        int LWORK = -1;
        int INFO  = -1;

        {// find optimal work-space
          execute(JOBZ, RANGE, UPLO, N_A, &X(0,0), LDA, VL, VU, IL, IU, &lambda[0], &V(0,0), LDZ, &WORK[0], LWORK, &IWORK[0], &IFAIL[0], INFO);

          LWORK = WORK[0];

          WORK.resize(LWORK);
        }

        M = execute(JOBZ, RANGE, UPLO, N_A, &X(0,0), LDA, VL, VU, IL, IU, &lambda[0], &V(0,0), LDZ, &WORK[0], LWORK, &IWORK[0], &IFAIL[0], INFO);
      }
    else
      {
        throw std::logic_error("divide and conquer not yet implemented");
      }

    return M;
  }

  template<typename scalartype>
  int GEEV<CPU>::execute(char JOBZ, char RANGE, char UPLO,
                         scalartype VL,
                         scalartype VU,
                         int IL,
                         int IU,
                         matrix<std::complex<scalartype>, CPU>& A,
                         dca::linalg::Vector<scalartype, CPU>& lambda,
                         matrix<std::complex<scalartype>, CPU>& V,
                         std::string                            alg_type)
  {
    if( JOBZ != 'N' and JOBZ != 'V')
      throw std::logic_error(__FUNCTION__);

    if( UPLO != 'U' and UPLO != 'L')
      throw std::logic_error(__FUNCTION__);

    matrix<scalartype, CPU> X(A);

    V = A;

    int N_A = A.size().first;

    int LDA = X.leadingDimension();
    int LDZ = V.leadingDimension();

    int M=-1;

    if(alg_type == "DEFAULT")
      {
        dca::linalg::Vector<std::complex<scalartype>, CPU> WORK(1);
        dca::linalg::Vector<scalartype, CPU> RWORK(7*N_A);

        dca::linalg::Vector<int, CPU> IWORK(5*N_A);
        dca::linalg::Vector<int, CPU> IFAIL(N_A);

        int LWORK  = -1;
        int LRWORK = -1;

        int INFO  = -1;

        {// find optimal work-space
          execute(JOBZ, RANGE, UPLO, N_A, &X(0,0), LDA, VL, VU, IL, IU, &lambda[0], &V(0,0), LDZ, &WORK[0], LWORK, &RWORK[0], LRWORK, &IWORK[0], &IFAIL[0], INFO);

          LWORK = WORK[0];

          WORK.resize(LWORK);
        }

        M = execute(JOBZ, RANGE, UPLO, N_A, &X(0,0), LDA, VL, VU, IL, IU, &lambda[0], &V(0,0), LDZ, &WORK[0], LWORK, &RWORK[0], LRWORK, &IWORK[0], &IFAIL[0], INFO);
      }
    else
      {
        throw std::logic_error("divide and conquer not yet implemented");
      }

    return M;
  }


  /**************************
   ***
   ***   GEEV LAPACK-routines
   ***
   **************************/

  void GEEV<CPU>::execute(char JOBVL, char JOBVR, int N, float* A,  int LDA, float* WR, float* WI, float* VL, int LDVL, float* VR, int LDVR, float* WORK, int LWORK,int INFO )
  {
    LAPACK::sgeev_(&JOBVL, &JOBVR, &N, A, &LDA, WR, WI, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  void GEEV<CPU>::execute(char JOBVL, char JOBVR, int N, double* A,  int LDA, double* WR, double* WI, double* VL, int LDVL, double* VR, int LDVR, double* WORK, int LWORK, int INFO )
  {
    LAPACK::dgeev_(&JOBVL, &JOBVR, &N, A, &LDA, WR, WI, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  void GEEV<CPU>::execute(char JOBVL, char JOBVR, int N, std::complex<float>* A, int LDA, std::complex<float>* W,
                          std::complex<float>* VL, int LDVL, std::complex<float>* VR, int LDVR,
                          std::complex<float>* WORK, int LWORK, float* RWORK, int INFO)
  {
    LAPACK::cgeev_(&JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR, &LDVR, WORK, &LWORK, RWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  void GEEV<CPU>::execute(char JOBVL, char JOBVR, int N, std::complex<double>* A, int LDA, std::complex<double>* W,
                          std::complex<double>* VL, int LDVL, std::complex<double>* VR, int LDVR,
                          std::complex<double>* WORK, int LWORK, double* RWORK, int INFO)
  {
    LAPACK::zgeev_(&JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR, &LDVR, WORK, &LWORK, RWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  /**************************
   ***
   ***   HEEV LAPACK-routines
   ***
   **************************/

  void GEEV<CPU>::execute(char JOBZ,  char UPLO,  int N, float*  A,  int LDA, float*  W, float* WORK, int LWORK)
  {
    int INFO = -1;

    LAPACK::ssyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  void GEEV<CPU>::execute(char JOBZ,  char UPLO,  int N, double* A,  int LDA, double* W, double* WORK,  int LWORK)
  {
    int INFO = -1;

    LAPACK::dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  void GEEV<CPU>::execute(char JOBZ, char UPLO, int N, std::complex<float>* A, int LDA, float*  W, std::complex<float>* WORK, int LWORK, float* RWORK)
  {
    int INFO = -1;

    LAPACK::cheev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, RWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  void GEEV<CPU>::execute(char JOBZ, char UPLO, int N, std::complex<double>* A, int LDA, double* W, std::complex<double>* WORK,  int LWORK, double* RWORK)
  {
    int INFO = -1;

    LAPACK::zheev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, RWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  /**************************
   ***
   ***   HEEVD LAPACK-routines
   ***
   **************************/

  void GEEV<CPU>::execute(char JOBZ,  char UPLO,  int N, float*  A,  int LDA, float*  W,
                          float*  WORK, int LWORK, int* IWORK, int LIWORK)
  {
    int INFO = -1;

    LAPACK::ssyevd_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, IWORK, &LIWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  void GEEV<CPU>::execute(char JOBZ,  char UPLO,  int N, double* A,  int LDA, double* W,
                          double* WORK, int LWORK, int* IWORK, int LIWORK)
  {
    int INFO = -1;

    LAPACK::dsyevd_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, IWORK, &LIWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  void GEEV<CPU>::execute(char JOBZ,  char UPLO, int N, std::complex<float>*  A,  int LDA, float*  W,
                          std::complex<float>* WORK , int LWORK,
                          float*               RWORK, int LRWORK,
                          int*                 IWORK, int LIWORK)
  {
    int INFO = -1;

    LAPACK::cheevd_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  void GEEV<CPU>::execute(char JOBZ,  char UPLO, int N, std::complex<double>* A,  int LDA, double* W,
                          std::complex<double>* WORK , int LWORK,
                          double*               RWORK, int LRWORK,
                          int*                  IWORK, int LIWORK)
  {
    int INFO = -1;

    LAPACK::zheevd_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  /****************************
   ***
   ***   HEEVX LAPACK-routines
   ***
   ****************************/

  int GEEV<CPU>::execute(char JOBZ,  char RANGE, char UPLO,  int N, float * A, int LDA,
                         float VL, float VU, int IL, int UL,
                         float* W, float* Z, int LDZ,
                         float* WORK, int LWORK, int* IWORK, int* IFAIL, int INFO )
  {
    int M = -1;

    char  tmp    = 'S';
    float ABSTOL = 2*LAPACK::slamch_(&tmp);

    LAPACK::ssyevx_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &UL,
                    &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, IWORK, IFAIL, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);

    return M;
  }

  int GEEV<CPU>::execute(char JOBZ,  char RANGE, char UPLO,  int N, double* A, int LDA,
                         double VL, double VU, int IL, int UL,
                         double* W, double* Z, int LDZ,
                         double* WORK, int LWORK, int* IWORK, int* IFAIL, int INFO )
  {
    int M = -1;

    char   tmp    = 'S';
    double ABSTOL = 2*LAPACK::dlamch_(&tmp);

    LAPACK::dsyevx_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &UL,
                    &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, IWORK, IFAIL, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);

    return M;
  }

  int GEEV<CPU>::execute(char JOBZ,  char RANGE, char UPLO, int N, std::complex<float>* A, int LDA,
                         float VL, float VU, int IL, int UL,
                         float* W, std::complex<float>* Z, int LDZ,
                         std::complex<float>* WORK, int LWORK, float* RWORK, int LRWORK, int* IWORK, int* IFAIL, int INFO )
  {
    int M = -1;

    char  tmp    = 'S';
    float ABSTOL = 2*LAPACK::slamch_(&tmp);

    LAPACK::cheevx_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &UL,
                    &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, RWORK, &LRWORK, IWORK, IFAIL, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);

    return M;
  }

  int GEEV<CPU>::execute(char JOBZ,  char RANGE, char UPLO, int N, std::complex<double>* A, int LDA,
                         double VL, double VU, int IL, int UL,
                         double* W, std::complex<double>* Z, int LDZ,
                         std::complex<double>* WORK, int LWORK, double* RWORK, int LRWORK, int* IWORK, int* IFAIL, int INFO )
  {
    int M = -1;

    char  tmp    = 'S';
    double ABSTOL = 2*LAPACK::dlamch_(&tmp);

    LAPACK::zheevx_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &UL,
                    &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, RWORK, &LRWORK, IWORK, IFAIL, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);

    return M;
  }


  /****************************
   ***
   ***   HEEVR LAPACK-routines
   ***
   ****************************/

  void GEEV<CPU>::execute(char JOBZ,  char RANGE, char UPLO,  int N, float* A, int LDA,
                          float VL, float VU, int IL, int UL, int M,
                          float* W, float* Z, int LDZ, int* ISUPPZ,
                          float* WORK,  int LWORK, int* IWORK, int LIWORK, int INFO)
  {
    char  tmp    = 'S';
    float ABSTOL = 2*LAPACK::slamch_(&tmp);

    LAPACK::ssyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &UL,
                    &ABSTOL, &M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }


  void GEEV<CPU>::execute(char JOBZ,  char RANGE, char UPLO,  int N, double* A, int LDA,
                          double VL, double VU, int IL, int UL, int M,
                          double* W, double* Z, int LDZ, int* ISUPPZ,
                          double* WORK, int LWORK, int* IWORK, int LIWORK, int INFO)
  {
    char  tmp    = 'S';
    double ABSTOL = 2*LAPACK::dlamch_(&tmp);

    LAPACK::dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &UL,
                    &ABSTOL, &M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }


  void GEEV<CPU>::execute(char JOBZ,  char RANGE, char UPLO,  int N, std::complex<float>* A, int LDA,
                          float VL, float VU, int IL, int UL, int M,
                          float* W, std::complex<float>* Z, int LDZ, int* ISUPPZ,
                          std::complex<float>* WORK, int LWORK, float* RWORK, int LRWORK, int* IWORK, int LIWORK, int INFO)
  {
    char  tmp    = 'S';
    float ABSTOL = 2*LAPACK::slamch_(&tmp);

    LAPACK::cheevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &UL,
                    &ABSTOL, &M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

  void GEEV<CPU>::execute(char JOBZ,  char RANGE, char UPLO,  int N, std::complex<double>* A, int LDA,
                          double VL, double VU, int IL, int UL, int M,
                          double* W, std::complex<double>* Z, int LDZ, int* ISUPPZ,
                          std::complex<double>* WORK, int LWORK, double* RWORK, int LRWORK, int* IWORK, int LIWORK, int INFO)
  {
    char  tmp    = 'S';
    double ABSTOL = 2*LAPACK::dlamch_(&tmp);

    LAPACK::zheevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &UL,
                    &ABSTOL, &M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO);

    if(INFO != 0)
      throw std::logic_error(__FUNCTION__);
  }

}

#endif
