//-*-C++-*-                                                                                                                                                                                                                                                                                                                                                        

#ifndef LINALG_GEEV_GPU_CU_H
#define LINALG_GEEV_GPU_CU_H

namespace LIN_ALG 
{
  namespace GPU_KERNELS_GEEV
  {
     /**************************
      ***   GEEV-routines
      **************************/

     void sgeev(char JOBVL,  char JOBVR, int N, float* A,  int LDA, float* WR, float* WI, float* VL, int LDVL, float* VR, int LDVR, float* WORK, int LWORK)
     {
       int INFO = -1;
       
       magma_sgeev(magma_vec_const(JOBVL), magma_vec_const(JOBVR), N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, &INFO);

       if(INFO != 0)
	 throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
     }

     void dgeev(char JOBVL,  char JOBVR, int N, double* A,  int LDA, double* WR, double* WI, double* VL, int LDVL, double* VR, int LDVR, double* WORK, int LWORK)
     {
       int INFO = -1;
       
       magma_dgeev(magma_vec_const(JOBVL), magma_vec_const(JOBVR), N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, &INFO);

       if(INFO != 0)
	 throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
     }

     void cgeev(char JOBVL,  char JOBVR, int N, std::complex<float>* A,  int LDA, std::complex<float>* W, std::complex<float>* VL, int LDVL, std::complex<float>* VR, int LDVR, std::complex<float>* WORK, int LWORK, float* RWORK)
     {
       int INFO = -1;
       
       magma_cgeev(magma_vec_const(JOBVL), magma_vec_const(JOBVR), N, (magmaFloatComplex*) A, LDA, (magmaFloatComplex*) W, (magmaFloatComplex*) VL, LDVL, (magmaFloatComplex*) VR, LDVR, (magmaFloatComplex*) WORK, LWORK, RWORK, &INFO);

       if(INFO != 0)
	 throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
     }

     void zgeev(char JOBVL,  char JOBVR, int N, std::complex<double>* A,  int LDA, std::complex<double>* W, std::complex<double>* VL, int LDVL, std::complex<double>* VR, int LDVR, std::complex<double>* WORK, int LWORK, double* RWORK)
     {
       int INFO = -1;
       
       magma_zgeev(magma_vec_const(JOBVL), magma_vec_const(JOBVR), N, (magmaDoubleComplex*) A, LDA, (magmaDoubleComplex*) W, (magmaDoubleComplex*) VL, LDVL, (magmaDoubleComplex*) VR, LDVR, (magmaDoubleComplex*) WORK, LWORK, RWORK, &INFO);

       if(INFO != 0)
	 throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
     }

     /**************************
      ***   HEEVD-routines
      **************************/

     void ssyevd( char JOBZ,  char UPLO,  int N, float* A, int LDA, float* W , float* WORK , int LWORK, int* IWORK, int LIWORK)
     {
       int INFO = -1;
       
       magma_ssyevd(magma_vec_const(JOBZ), magma_uplo_const(UPLO), N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, &INFO);

       if(INFO != 0)
	 throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
     }

     void dsyevd( char JOBZ,  char UPLO,  int N, double* A, int LDA, double* W, double* WORK, int LWORK, int* IWORK, int LIWORK)
     {
       int INFO = -1;
       
       magma_dsyevd(magma_vec_const(JOBZ), magma_uplo_const(UPLO), N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, &INFO);

       if(INFO != 0)
	 throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
     }

    void cheevd(char JOBZ,  char UPLO, int N, std::complex<float>* A,  int LDA, float* W, std::complex<float>* WORK, int LWORK, float* RWORK, int LRWORK, int* IWORK, int LIWORK)
    {
      int INFO = -1;
      
      magma_cheevd(magma_vec_const(JOBZ), magma_uplo_const(UPLO), N, (magmaFloatComplex*) A, LDA, W, (magmaFloatComplex*) WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, &INFO);
      
      if(INFO != 0)
	throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }
    
    void zheevd(char JOBZ, char UPLO, int N, std::complex<double>* A, int LDA, double* W, std::complex<double>* WORK, int LWORK, double* RWORK, int LRWORK, int* IWORK, int LIWORK)
    {
      int INFO = -1;
      
      magma_zheevd(magma_vec_const(JOBZ), magma_uplo_const(UPLO), N, (magmaDoubleComplex*) A, LDA, W, (magmaDoubleComplex*) WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, &INFO);
      
      if(INFO != 0)
	throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif
    }

     /**************************
      ***   HEEVR-routines
      **************************/

     int ssyevdx(char JOBZ,  char RANGE, char UPLO,  int N, float * A, int LDA, 
		 float VL, float VU, int IL, int UL,
		 float* W, float* WORK,  int LWORK, int* IWORK, int LIWORK)
     {
       int M        =  0;
       int INFO     = -1;

       magma_ssyevdx(magma_vec_const(JOBZ), magma_range_const(RANGE), magma_uplo_const(UPLO), N, A, LDA, VL, VU, IL, UL, &M, W, WORK, LWORK, IWORK, LIWORK, &INFO);

       if(INFO != 0)
	 throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif

       return M;
     }

     int dsyevdx(char JOBZ,  char RANGE, char UPLO,  int N, double* A, int LDA, 
		 double VL, double VU, int IL, int UL,
		 double* W, double* WORK,  int LWORK, int* IWORK, int LIWORK)
     {
       int M         =  0;
       int INFO      = -1;
       
       magma_dsyevdx(magma_vec_const(JOBZ), magma_range_const(RANGE), magma_uplo_const(UPLO), N, A, LDA, VL, VU, IL, UL, &M, W, WORK, LWORK, IWORK, LIWORK, &INFO);

       if(INFO != 0)
	 throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif

       return M;
     }

     int cheevdx(char JOBZ,  char RANGE, char UPLO,  int N, std::complex<float>* A, int LDA, 
		 float VL, float VU, int IL, int UL,
		 float* W, std::complex<float>* WORK, int LWORK, float* RWORK,  int LRWORK, int* IWORK, int LIWORK)
     {
       int M      =  0;
       int INFO = -1;

       magma_cheevdx(magma_vec_const(JOBZ), magma_range_const(RANGE), magma_uplo_const(UPLO), N, (magmaFloatComplex*) A, LDA, VL, VU, IL, UL, &M, W, (magmaFloatComplex*) WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, &INFO);

       if(INFO != 0)
	 throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif

       return M;
     }

     int zheevdx(char JOBZ,  char RANGE, char UPLO,  int N, std::complex<double>* A, int LDA, 
		 double VL, double VU, int IL, int UL,
		 double* W, std::complex<double>* WORK,  int LWORK, double* RWORK, int LRWORK, int* IWORK, int LIWORK)
     {
       int M    = 0;
       int INFO = -1;

       magma_zheevdx(magma_vec_const(JOBZ), magma_range_const(RANGE), magma_uplo_const(UPLO), N, (magmaDoubleComplex*) A, LDA, VL, VU, IL, UL, &M, W, (magmaDoubleComplex*) WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, &INFO);

       if(INFO != 0)
	 throw std::logic_error(__FUNCTION__);

#ifdef DEBUG_CUDA
       cuda_check_for_errors(__FUNCTION__, __FILE__, __LINE__);
#endif

       return M;
     }
  }

}


#endif
