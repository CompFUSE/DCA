// //-*-C++-*-

// #ifndef ANY_GETRS_PLAN_H
// #define ANY_GETRS_PLAN_H

// /*!
//  *   \author Peter Staar
//  */
// template<typename scalartype>
// class getrs_plan<scalartype, ANY_LIBRARY>
// {
// public:
    
//   getrs_plan(int n, int nrhs, int lda, int ldb);
//   ~getrs_plan();
    
//   inline void execute_plan(LINEAR_ALGEBRA_LIBRARY_TYPE LAL_t);
    
// private:

//   template<typename whatever_t>
//   inline void set_values(whatever_t& whatever_ref);
    
// public:
    
//   char TRANS;
//   int  N;
//   int  NRHS;
//   int  LDA;
//   int  LDB;
//   int  INFO;

//   scalartype*  Matrix_A;
//   scalartype*  Matrix_B;
//   int*             IPIV;

// };

// template<typename scalartype>
// getrs_plan<scalartype, ANY_LIBRARY>::getrs_plan(int n, int nrhs, int lda, int ldb):
//   TRANS('N'),
//   N(n),
//   NRHS(nrhs),
//   LDA(lda),
//   LDB(ldb)
// {
//   IPIV = new int[N];
  
//   for(int i=0; i<N; i++)
//     IPIV[i] = i+1;
// }

// template<typename scalartype>
// getrs_plan<scalartype, ANY_LIBRARY>::~getrs_plan()
// {
//   delete [] IPIV;
// }

// template<typename scalartype>
// template<typename whatever_t>
// void getrs_plan<scalartype, ANY_LIBRARY>::set_values(whatever_t& whatever_ref)
// {
//   whatever_ref.TRANS = TRANS;

//   whatever_ref.N    = N;
//   whatever_ref.NRHS = NRHS;

//   whatever_ref.LDA = LDA;
//   whatever_ref.LDB = LDB;

//   whatever_ref.INFO = INFO;

//   whatever_ref.Matrix_A = &Matrix_A[0];
//   whatever_ref.Matrix_B = &Matrix_B[0];

//   memcpy(&whatever_ref.IPIV[0], &IPIV[0], sizeof(int)*N);
// }

// template<typename scalartype>
// void getrs_plan<scalartype, ANY_LIBRARY>::execute_plan(LINEAR_ALGEBRA_LIBRARY_TYPE LAL_t)
// {
//   switch(LAL_t)
//     {
//     case LAPACK_LIBRARY:
//       {
// 	getrs_plan<scalartype, LAPACK_LIBRARY> lapack_getrs_plan(N, NRHS, LDA, LDB);
// 	set_values(lapack_getrs_plan);
// 	lapack_getrs_plan.execute_plan();
//       }
//       break;

// #ifdef USE_MAGMA_LIBRARY
//     case MAGMA_LIBRARY:
//       {
// 	getrs_plan<scalartype, MAGMA_LIBRARY> magma_getrs_plan();
// 	set_values(magma_getrs_plan);
// 	magma_getrs_plan.execute_plan();
//       }
//       break;
// #endif

//     default:
//       throw std::logic_error(__FUNCTION__);
//     }
// }


// #endif
