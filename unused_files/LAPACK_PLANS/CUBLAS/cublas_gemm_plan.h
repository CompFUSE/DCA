// //-*-C++-*-

// #ifndef CUBLAS_MATRIX_MATRIX_PLAN_H
// #define CUBLAS_MATRIX_MATRIX_PLAN_H

// /*!
//  *   \author R. Solca, P. Staar
//  */
// template<typename scalartype>
// class gemm_plan<scalartype, CUBLAS_LIBRARY>
// {
//   typedef typename GET_CUBLAS_TYPE<scalartype>::type cuda_scalartype;

// public:
    
//   gemm_plan();
    
//   gemm_plan(int n);
    
//   gemm_plan(int m, int k, int n);
    
//   ~gemm_plan();
    
//   void  execute_plan();
    
// private:

//   void allocate_gpu();
 
//   void deallocate_gpu();
    
// public:
    
//   char TRANSA;
//   char TRANSB;
    
//   int M;
//   int N;
//   int K;
    
//   int LDA;
//   int LDB;
//   int LDC;
    
//   scalartype alpha;
//   scalartype beta;
    
//   scalartype* A;
//   scalartype* B;
//   scalartype* C;
    
//   cuda_scalartype* DA;
//   cuda_scalartype* DB;
//   cuda_scalartype* DC;
// };

// template<typename scalartype>
// gemm_plan<scalartype, CUBLAS_LIBRARY>::gemm_plan():
//   TRANSA('N'),
//   TRANSB('N'),

//   alpha(1.),
//   beta(0.)
// {}

// template<typename scalartype>
// gemm_plan<scalartype, CUBLAS_LIBRARY>::gemm_plan(int n):
//   TRANSA('N'),
//   TRANSB('N'),

//   M(n),
//   N(n),
//   K(n),

//   LDA(n),
//   LDB(n),
//   LDC(n),

//   alpha(1.),
//   beta(0.)
// {}

// template<typename scalartype>
// gemm_plan<scalartype, CUBLAS_LIBRARY>::gemm_plan(int m, int k, int n):
//   TRANSA('N'),
//   TRANSB('N'),

//   M(m),
//   N(n),
//   K(k),

//   LDA(m),
//   LDB(k),
//   LDC(m),

//   alpha(1.),
//   beta(0.)
// {}

// template<typename scalartype>
// gemm_plan<scalartype, CUBLAS_LIBRARY>::~gemm_plan()
// {}

// template<typename scalartype>
// void gemm_plan<scalartype, CUBLAS_LIBRARY>::allocate_gpu()
// {
//   cudaError_t ret[3];
  
//   ret[0] = cudaMalloc( (void**)&DA, K*LDA*sizeof(typename GET_CUBLAS_TYPE<scalartype>::type) );
//   ret[1] = cudaMalloc( (void**)&DB, N*LDB*sizeof(typename GET_CUBLAS_TYPE<scalartype>::type) );
//   ret[2] = cudaMalloc( (void**)&DC, N*LDC*sizeof(typename GET_CUBLAS_TYPE<scalartype>::type) );
    
//   if( ret[0] != cudaSuccess || ret[1] != cudaSuccess || ret[2] != cudaSuccess){
//     std::cout << "NOT ENOUGH GPU MEMORY " << (K*LDA+N*LDB+N*LDC)*sizeof(typename GET_CUBLAS_TYPE<scalartype>::type) 
// 	      << " ret code: " << (ret[0] != cudaSuccess) << (ret[1] != cudaSuccess) << (ret[2] != cudaSuccess)
// 	      << std::endl;
//     abort();
//   }
// }

// template<typename scalartype>
// void gemm_plan<scalartype, CUBLAS_LIBRARY>::deallocate_gpu()
// {
//   cudaFree(DA);
//   cudaFree(DB);
//   cudaFree(DC);
// }

// template<>
// void gemm_plan<float, CUBLAS_LIBRARY>::execute_plan()
// {
//   allocate_gpu();
    
//   cublasSetMatrix(M, K, sizeof(cuda_scalartype), A, LDA, DA, LDA);
//   cublasSetMatrix(K, N, sizeof(cuda_scalartype), B, LDB, DB, LDB);

//   if (beta != 0.)
//     cublasSetMatrix(M, N, sizeof(cuda_scalartype), C, LDC, DC, LDC);
    
//   cublasSgemm(TRANSA, TRANSB, M, N, K, alpha, DA, LDA, DB, LDB, beta, DC, LDC);
    
//   cublasGetMatrix(M, N, sizeof(cuda_scalartype), DC, LDC, C, LDC);
    
//   deallocate_gpu();
// }

// template<>
// void gemm_plan<double, CUBLAS_LIBRARY>::execute_plan()
// {
//   allocate_gpu();
  
//   cublasSetMatrix(M, K, sizeof(cuda_scalartype), A, LDA, DA, LDA);
//   cublasSetMatrix(K, N, sizeof(cuda_scalartype), B, LDB, DB, LDB);
  
//   if (beta != 0.)
//     cublasSetMatrix(M, N, sizeof(cuda_scalartype), C, LDC, DC, LDC);
  
//   cublasDgemm(TRANSA, TRANSB, M, N, K, alpha, DA, LDA, DB, LDB, beta, DC, LDC);
  
//   cublasGetMatrix(M, N, sizeof(cuda_scalartype), DC, LDC, C, LDC);
  
//   deallocate_gpu();
// }

// template<>
// void gemm_plan<std::complex<float>, CUBLAS_LIBRARY>::execute_plan()
// {
//   allocate_gpu();
  
//   cublasSetMatrix(M, K, sizeof(cuda_scalartype), A, LDA, DA, LDA);
//   cublasSetMatrix(K, N, sizeof(cuda_scalartype), B, LDB, DB, LDB);

//   cuda_scalartype alpha_;
//   alpha_.x = alpha.real();
//   alpha_.y = alpha.imag();

//   cuda_scalartype beta_;
//   beta_.x = beta.real();
//   beta_.y = beta.imag();

//   if (beta.real()==0 && beta.imag()==0)
//     cublasSetMatrix(M, N, sizeof(cuda_scalartype), C, LDC, DC, LDC);
  
//   cublasCgemm(TRANSA, TRANSB, M, N, K, alpha_, DA, LDA, DB, LDB, beta_, DC, LDC);

//   cublasGetMatrix(M, N, sizeof(cuda_scalartype), DC, LDC, C, LDC);
    
//   deallocate_gpu();
// }

// template<>
// void gemm_plan<std::complex<double>, CUBLAS_LIBRARY>::execute_plan()
// {
//   allocate_gpu();
  
//   cublasSetMatrix(M, K, sizeof(cuda_scalartype), A, LDA, DA, LDA);
//   cublasSetMatrix(K, N, sizeof(cuda_scalartype), B, LDB, DB, LDB);

//   cuda_scalartype alpha_;
//   alpha_.x = alpha.real();
//   alpha_.y = alpha.imag();
  
//   cuda_scalartype beta_;
//   beta_.x = beta.real();
//   beta_.y = beta.imag();
  
//   if (beta.real()==0 && beta.imag()==0)
//     cublasSetMatrix(M, N, sizeof(cuda_scalartype), C, LDC, DC, LDC);

//   cublasZgemm(TRANSA, TRANSB, M, N, K, alpha_, DA, LDA, DB, LDB, beta_, DC, LDC);

//   cublasGetMatrix(M, N, sizeof(cuda_scalartype), DC, LDC, C, LDC);
  
//   deallocate_gpu();
// }

// #endif
