//-*-C++-*-

#ifndef CUBLAS_GEMM_H
#define CUBLAS_GEMM_H

namespace CUDA_GPU
{

  /*!
   *   \author R. Solca, P. Staar
   */
  
  template<typename scalartype>
  void cugemm_t(char& TRANSA, char& TRANSB, int& M, int& N, int& K, scalartype& alpha,
		scalartype* DA, int& LDDA, scalartype* DB, int& LDDB, scalartype& beta,
		scalartype* DC, int& LDDC)
  {
    cout << __PRETTY_FUNCTION__ << endl;
    throw std::logic_error(__FUNCTION__);
  }
  
  template<>
  void cugemm_t<float>(char& TRANSA, char& TRANSB, int& M, int& N, int& K, float& alpha,
		       float* DA, int& LDDA, float* DB, int& LDDB, float& beta,
		       float* DC, int& LDDC)
  {
    cublasSgemm(TRANSA, TRANSB, M, N, K, alpha, DA, LDDA, DB, LDDB, beta, DC, LDDC);
  }

  template<>
  void cugemm_t<double>(char& TRANSA, char& TRANSB, int& M, int& N, int& K, double& alpha,
			double* DA, int& LDDA, double* DB, int& LDDB, double& beta,
			double* DC, int& LDDC)
  {
    cublasDgemm(TRANSA, TRANSB, M, N, K, alpha, DA, LDDA, DB, LDDB, beta, DC, LDDC);
  }

  template<>
  void cugemm_t<cuComplex>(char& TRANSA, char& TRANSB, int& M, int& N, int& K, cuComplex& alpha,
			   cuComplex* DA, int& LDDA, cuComplex* DB, int& LDDB, cuComplex& beta,
			   cuComplex* DC, int& LDDC)
  {
    cublasCgemm(TRANSA, TRANSB, M, N, K, alpha, DA, LDDA, DB, LDDB, beta, DC, LDDC);
  }
  
  template<>
  void cugemm_t<cuDoubleComplex>(char& TRANSA, char& TRANSB, int& M, int& N, int& K, cuDoubleComplex& alpha,
				 cuDoubleComplex* DA, int& LDDA, cuDoubleComplex* DB, int& LDDB, cuDoubleComplex& beta,
				 cuDoubleComplex* DC, int& LDDC)
  {
    cublasZgemm(TRANSA, TRANSB, M, N, K, alpha, DA, LDDA, DB, LDDB, beta, DC, LDDC);
  }
  
  template<typename scalartype>
  void cublas_gemm_t(char& TRANSA, char& TRANSB, int& M, int& N, int& K, scalartype& alpha,
		     scalartype* A, int& LDA, scalartype* B, int& LDB, scalartype& beta,
		     scalartype* C, int& LDC)
  {
    //std::cout << "cublas_gemm " << M << " " << N << " " << K << std::endl;
    
    typedef typename GET_CUBLAS_TYPE<scalartype>::type cuda_scalartype;
    
    cuda_scalartype *DA = NULL;
    cuda_scalartype *DB = NULL;
    cuda_scalartype *DC = NULL;
    
    int LDDA=LDA;
    int LDDB=LDB;
    int LDDC=LDC;

    allocate_gpu(&DA, LDDA, K);
    allocate_gpu(&DB, LDDB, N);
    allocate_gpu(&DC, LDDC, N);
    
    cublasSetMatrix(M, K, sizeof(cuda_scalartype), A, LDA, DA, LDDA);
    cublasSetMatrix(K, N, sizeof(cuda_scalartype), B, LDB, DB, LDDB);
    
    if (beta != scalartype(0))
        cublasSetMatrix(M, N, sizeof(cuda_scalartype), C, LDC, DC, LDDC);
    
    cuda_scalartype alpha_;
    CAST_CUBLAS_TYPE<scalartype, cuda_scalartype>::execute(alpha, alpha_);
    cuda_scalartype beta_;
    CAST_CUBLAS_TYPE<scalartype, cuda_scalartype>::execute(beta, beta_);
    
    cugemm_t(TRANSA, TRANSB, M, N, K, alpha_, DA, LDDA, DB, LDDB, beta_, DC, LDDC);
    
    //std::cout << "cugemm error: " << cudaDeviceSynchronize() << std::endl;
    
    cublasGetMatrix(M, N, sizeof(cuda_scalartype), DC, LDDC, C, LDC);
    
    deallocate_gpu(DA);
    deallocate_gpu(DB);
    deallocate_gpu(DC);
  }
}
#endif
