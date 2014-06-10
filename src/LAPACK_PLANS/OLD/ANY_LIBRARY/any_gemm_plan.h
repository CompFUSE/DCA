// //-*-C++-*-

// #ifndef ANY_MATRIX_MATRIX_PLAN_H
// #define ANY_MATRIX_MATRIX_PLAN_H

// /*!
//  *   \author P. Staar
//  */
// template<typename scalartype>
// class gemm_plan<scalartype, ANY_LIBRARY>
// {
// public:
    
//   gemm_plan();
    
//   gemm_plan(int n);
    
//   gemm_plan(int m, int k, int n);
    
//   ~gemm_plan();
    
//   inline void execute_plan(LINEAR_ALGEBRA_LIBRARY_TYPE LAL_t);
    
// private:

//   template<typename whatever_t>
//   inline void set_values(whatever_t& whatever_ref);
    
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
// };

// template<typename scalartype>
// gemm_plan<scalartype, ANY_LIBRARY>::gemm_plan():
//   TRANSA('N'),
//   TRANSB('N'),

//   alpha(1.),
//   beta(0.)
// {}

// template<typename scalartype>
// gemm_plan<scalartype, ANY_LIBRARY>::gemm_plan(int n):
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
// gemm_plan<scalartype, ANY_LIBRARY>::gemm_plan(int m, int k, int n):
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
// gemm_plan<scalartype, ANY_LIBRARY>::~gemm_plan()
// {}

// template<typename scalartype>
// template<typename whatever_t>
// void gemm_plan<scalartype, ANY_LIBRARY>::set_values(whatever_t& whatever_ref)
// {
//   whatever_ref.TRANSA = TRANSA;
//   whatever_ref.TRANSB = TRANSB;

//   whatever_ref.M = M;
//   whatever_ref.N = N;
//   whatever_ref.K = K;

//   whatever_ref.LDA = LDA;
//   whatever_ref.LDB = LDB;
//   whatever_ref.LDC = LDC;

//   whatever_ref.alpha = alpha;
//   whatever_ref.beta  = beta;

//   whatever_ref.A = &A[0];
//   whatever_ref.B = &B[0];
//   whatever_ref.C = &C[0];
// }

// template<typename scalartype>
// void gemm_plan<scalartype, ANY_LIBRARY>::execute_plan(LINEAR_ALGEBRA_LIBRARY_TYPE LAL_t)
// {
//   switch(LAL_t)
//     {
//     case BLAS_LIBRARY:
//       {
// 	gemm_plan<scalartype, BLAS_LIBRARY> blas_gemm_plan;
// 	set_values(blas_gemm_plan);
// 	blas_gemm_plan.execute_plan();
//       }
//       break;

// #ifdef USE_CUBLAS_LIBRARY
//     case CUBLAS_LIBRARY:
//       {
// 	gemm_plan<scalartype, CUBLAS_LIBRARY> cublas_gemm_plan();
// 	set_values(cublas_gemm_plan);
// 	cublas_gemm_plan.execute_plan();
//       }
//       break;
// #endif

// #ifdef USE_MAGMA_LIBRARY
//     case MAGMA_LIBRARY:
//       {
// 	gemm_plan<scalartype, MAGMA_LIBRARY> magma_gemm_plan();
// 	set_values(magma_gemm_plan);
// 	magma_gemm_plan.execute_plan();
//       }
//       break;
// #endif

//     default:
//       throw std::logic_error(__FUNCTION__);
//     }
// }


// #endif
