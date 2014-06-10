// //-*-C++-*-

// #ifndef MATRIX_MATRIX_PLAN_H
// #define MATRIX_MATRIX_PLAN_H

// /*!
//  * matrix_matrix_multiplication_plan.h
//  *
//  *      Author: peter staar
//  */
// template<typename scalartype, LINEAR_ALGEBRA_LIBRARY_TYPE linear_algebra_library=BLAS_LIBRARY>
// class gemm_plan
// {

//  public:

//   gemm_plan();

//   gemm_plan(int n);

//   gemm_plan(int m, int k, int n);

//   ~gemm_plan();

//   void  execute_plan();

//  public:

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

// // *  Purpose
// // *  =======
// // *
// // *  ZGEMM  performs one of the matrix-matrix operations
// // *
// // *     C := alpha*op( A )*op( B ) + beta*C,
// // *
// // *  where  op( X ) is one of
// // *
// // *     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
// // *
// // *  alpha and beta are scalars, and A, B and C are matrices, with op( A )
// // *  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
// // *
// // *  Arguments
// // *  ==========
// // *
// // *  TRANSA - CHARACTER*1.
// // *           On entry, TRANSA specifies the form of op( A ) to be used in
// // *           the matrix multiplication as follows:
// // *
// // *              TRANSA = 'N' or 'n',  op( A ) = A.
// // *
// // *              TRANSA = 'T' or 't',  op( A ) = A'.
// // *
// // *              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
// // *
// // *           Unchanged on exit.
// // *
// // *  TRANSB - CHARACTER*1.
// // *           On entry, TRANSB specifies the form of op( B ) to be used in
// // *           the matrix multiplication as follows:
// // *
// // *              TRANSB = 'N' or 'n',  op( B ) = B.
// // *
// // *              TRANSB = 'T' or 't',  op( B ) = B'.
// // *
// // *              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
// // *
// // *           Unchanged on exit.
// // *
// // *  M      - INTEGER.
// // *           On entry,  M  specifies  the number  of rows  of the  matrix
// // *           op( A )  and of the  matrix  C.  M  must  be at least  zero.
// // *           Unchanged on exit.
// // *
// // *  N      - INTEGER.
// // *           On entry,  N  specifies the number  of columns of the matrix
// // *           op( B ) and the number of columns of the matrix C. N must be
// // *           at least zero.
// // *           Unchanged on exit.
// // *
// // *  K      - INTEGER.
// // *           On entry,  K  specifies  the number of columns of the matrix
// // *           op( A ) and the number of rows of the matrix op( B ). K must
// // *           be at least  zero.
// // *           Unchanged on exit.

// template<typename scalartype, LINEAR_ALGEBRA_LIBRARY_TYPE linear_algebra_library>
// gemm_plan<scalartype, linear_algebra_library>::gemm_plan():
//   TRANSA('N'),
//   TRANSB('N'),

//   alpha(1.),
//   beta(0.)
// {}

// template<typename scalartype, LINEAR_ALGEBRA_LIBRARY_TYPE linear_algebra_library>
// gemm_plan<scalartype, linear_algebra_library>::gemm_plan(int n):
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

// template<typename scalartype, LINEAR_ALGEBRA_LIBRARY_TYPE linear_algebra_library>
// gemm_plan<scalartype, linear_algebra_library>::gemm_plan(int m, int k, int n):
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

// template<typename scalartype, LINEAR_ALGEBRA_LIBRARY_TYPE linear_algebra_library>
// gemm_plan<scalartype, linear_algebra_library>::~gemm_plan()
// {}

// template<typename scalartype, LINEAR_ALGEBRA_LIBRARY_TYPE linear_algebra_library>
// void gemm_plan<scalartype, linear_algebra_library>::execute_plan()
// {
//   cout << __PRETTY_FUNCTION__ << endl;
//   throw std::logic_error(__FUNCTION__);
// }

// template<>
// void gemm_plan<float, BLAS_LIBRARY>::execute_plan()
// {
//   BLAS::sgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
// }

// template<>
// void gemm_plan<double, BLAS_LIBRARY>::execute_plan()
// {
//   BLAS::dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
// }

// template<>
// void gemm_plan<std::complex<float>, BLAS_LIBRARY>::execute_plan()
// {
//   BLAS::cgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
// }

// template<>
// void gemm_plan<std::complex<double>, BLAS_LIBRARY>::execute_plan()
// {
//   BLAS::zgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
// }

// #endif
