//-*-C++-*-

/*
 * matrix_vector_multiplication_plan.h
 *
 *      Author: peter staar
 */

/*#include <cmath>
#include <cstdlib> 
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

using namespace std;*/


#ifndef BANDED_MATRIX_VECTOR_PLAN_H_
#define BANDED_MATRIX_VECTOR_PLAN_H_


template<typename scalartype>
class gbmv_plan
{

 public:

  gbmv_plan(int n);
  gbmv_plan(int m, int n);
  ~gbmv_plan();

  void execute_plan();

//   static void test_speed();

public:

  char TRANS;
  int M; // --> matrix-rows 
  int N; // --> matrix-columns 
  int KL; // number of sub-diagonals --> 0 for diagonal matrix
  int KU; // number of super-diagonals --> 0 for diagonal matrix
 
  scalartype   alpha;
  scalartype*  A; // -> matrix
  int LDA;

  scalartype* X; // -> source-vector
  int INCX;

  scalartype   beta;
  scalartype* Y; // -> target-vector
  int INCY;
};

// extern "C" void sgbmv_(const char* TRANS, const int* M, const int* N, const int* KL, const int* KU,
// 		      const float* alpha, const float* A, const int* LDA, 
// 		      const float* X, const int* INCX, 
// 		      const float* beta, float* Y, const int* INCY);

// extern "C" void dgbmv_(const char* TRANS, const int* M, const int* N, const int* KL, const int* KU, 
// 		      const double* alpha, const double* A, const int* LDA, 
// 		      const double* X, const int* INCX, 
// 		      const double* beta, double* Y, const int* INCY);

// extern "C" void cgbmv_(const char* TRANS, const int* M, const int* N, const int* KL, const int* KU,
// 		      const std::complex<float>* alpha, const std::complex<float>* A, const int* LDA,
// 		      const std::complex<float>* X, const int* INCX,
// 		      const std::complex<float>* beta, std::complex<float>* Y, const int* INCY);

// extern "C" void zgbmv_(const char* TRANS, const int* M, const int* N, const int* KL, const int* KU,
// 		      const std::complex<double>* alpha, const std::complex<double>* A, const int* LDA,
// 		      const std::complex<double>* X, const int* INCX, 
// 		      const std::complex<double>* beta, std::complex<double>* Y, const int* INCY);

template<typename scalartype>
gbmv_plan<scalartype>::gbmv_plan(int n):
  TRANS('N'),
  M(n),
  N(n),
  KL(0),
  KU(0),
  alpha(1.),
  LDA(1),
  INCX(1),
  beta(0.),
  INCY(1)
{}

template<typename scalartype>
gbmv_plan<scalartype>::gbmv_plan(int m, int n):
  TRANS('N'),
  M(m),
  N(n),
  KL(0),
  KU(0),
  alpha(1.),
  LDA(1),
  INCX(1),
  beta(0.),
  INCY(1)
{}

template<typename scalartype>
gbmv_plan<scalartype>::~gbmv_plan()
{}

template<typename scalartype>
void gbmv_plan<scalartype>::execute_plan()
{
  cout << __PRETTY_FUNCTION__ << endl;
  assert(false);
}

template<>
void gbmv_plan<float>::execute_plan()
{
  BLAS::sgbmv_(&TRANS, &M, &N, &KL, &KU, &alpha, A, &LDA, X, &INCX, &beta, Y, &INCY);
}

template<>
void gbmv_plan<double>::execute_plan()
{
  BLAS::dgbmv_(&TRANS, &M, &N, &KL, &KU, &alpha, A, &LDA, X, &INCX, &beta, Y, &INCY);
}

template<>
void gbmv_plan<std::complex<float> >::execute_plan()
{
  BLAS::cgbmv_(&TRANS, &M, &N, &KL, &KU, &alpha, A, &LDA, X, &INCX, &beta, Y, &INCY);
}

template<>
void gbmv_plan<std::complex<double> >::execute_plan()
{
  BLAS::zgbmv_(&TRANS, &M, &N, &KL, &KU, &alpha, A, &LDA, X, &INCX, &beta, Y, &INCY);
}

// template<typename scalartype>
// void gbmv_plan<scalartype>::test_speed()
// {
//   for(int i=100; i<2100; i += 100)
//     {
//       cout << i << "\t";

//       scalartype* m = new scalartype[i*i]; 
//       scalartype* v = new scalartype[i];
//       scalartype* M = new scalartype[i*i];

//       clock_t start_1 = clock();
//       gbmv_plan<scalartype> gbmv_pln(i);
//       for(int j=0; j<i; j++)
// 	{
// 	  gbmv_pln.A = v;
// 	  gbmv_pln.X = &m[i];
// 	  gbmv_pln.Y = &M[i];

// 	  gbmv_pln.execute_plan();
// 	}
//       clock_t end_1 = clock();

//       clock_t start_2 = clock();
//       for(int c=0; c<i; c++)
// 	for(int r=0; r<i; r++)
// 	  M[r + i*c] = v[r] * m[r + i*c];
//       clock_t end_2 = clock();

//       cout << double(int(end_1)-int(start_1))/double(CLOCKS_PER_SEC) << "\t" 
// 	   << double(int(end_2)-int(start_2))/double(CLOCKS_PER_SEC) << "\n";

//       delete [] m;
//       delete [] M;
//       delete [] v;
//     }
// }

#endif
