//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GESVD_GPU_CU_H
#define LINALG_GESVD_GPU_CU_H

namespace LIN_ALG {

  namespace GPU_KERNELS_GESVD
  {
    /**************************
     ***   GESVD-routines
     **************************/
    
    void sgesvd(char JOBU, char JOBVT, int M, int N, float*  A , int LDA, float*  S, float*  U , int LDU, float*  VT, int LDVT, float*  WORK, int LWORK)
    {
      int INFO = -1;

      magma_sgesvd(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, &INFO);

      if(INFO != 0)
	throw std::logic_error(__FUNCTION__);
    }
    
    void dgesvd(char JOBU, char JOBVT, int M, int N, double* A , int LDA, double* S, double* U , int LDU, double* VT, int LDVT, double* WORK, int LWORK)
    {
      int INFO = -1;

      magma_dgesvd(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, &INFO);

      if(INFO != 0)
	throw std::logic_error(__FUNCTION__);      
    }

    void cgesvd(char JOBU, char JOBVT, int M, int N, std::complex<float>*  A, int LDA, float*  S, std::complex<float>*  U, int LDU, std::complex<float>* VT, int LDVT, std::complex<float>*  WORK, int LWORK, float*  RWORK)
    {
      int INFO = -1;

      magma_cgesvd(JOBU, JOBVT, M, N, (magmaFloatComplex*) A, LDA, S, (magmaFloatComplex*) U, LDU, (magmaFloatComplex*) VT, LDVT, (magmaFloatComplex*) WORK, LWORK, RWORK, &INFO);
      
      if(INFO != 0)
	throw std::logic_error(__FUNCTION__);            
    }

    void zgesvd(char JOBU, char JOBVT, int M, int N, std::complex<double>* A, int LDA, double* S, std::complex<double>* U, int LDU, std::complex<double>* VT, int LDVT, std::complex<double>* WORK, int LWORK, double* RWORK)
    {
      int INFO = -1;

      magma_zgesvd(JOBU, JOBVT, M, N, (magmaDoubleComplex*) A, LDA, S, (magmaDoubleComplex*) U, LDU, (magmaDoubleComplex*) VT, LDVT, (magmaDoubleComplex*) WORK, LWORK, RWORK, &INFO);

      if(INFO != 0)
	throw std::logic_error(__FUNCTION__);                  
    }
    
    /**************************
     ***   GESVDD-routines
     **************************/
    
    /*
    void sgesdd(char JOBZ, int M, int N, float*  A, int LDA, float*  S, float*  U, int LDU, float*  VT, int LDVT, float*  WORK, int LWORK, int* IWORK)
    {
      int INFO = -1;

      magma_sgesdd(JOBZ,M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO);

      if(INFO != 0)
	throw std::logic_error(__FUNCTION__);                  
    }
    
    void dgesdd(char JOBZ, int M, int N, double* A, int LDA, double* S, double* U, int LDU, double* VT, int LDVT, double* WORK, int LWORK, int* IWORK)
    {
      int INFO = -1;

      magma_dgesdd(JOBZ,M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO);

      if(INFO != 0)
	throw std::logic_error(__FUNCTION__);                  
    }
    
    void cgesdd(char JOBZ, int M, int N, complex_float_t*  A, int LDA, float*  S, complex_float_t*  U, int LDU, complex_float_t*  VT, int LDVT, complex_float_t*  WORK, int LWORK, float*  RWORK, int* IWORK)
    {
      int INFO = -1;

      magma_cgesdd(JOBZ,M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, INFO);

      if(INFO != 0)
	throw std::logic_error(__FUNCTION__);                  
    }

    void zgesdd(char JOBZ, int M, int N, complex_double_t* A, int LDA, double* S, complex_double_t* U, int LDU, complex_double_t* VT, int LDVT, complex_double_t* WORK, int LWORK, double* RWORK, int* IWORK)
    {
      int INFO = -1;

      magma_zgesdd(JOBZ,M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, INFO);

      if(INFO != 0)
	throw std::logic_error(__FUNCTION__);                  
    }
    */
  }

}

#endif
