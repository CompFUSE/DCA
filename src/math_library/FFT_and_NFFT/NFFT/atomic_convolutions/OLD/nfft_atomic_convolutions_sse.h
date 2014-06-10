//-*-C++-*-

#ifndef MATH_ALGORITHMS_NFFT_ATOMIC_CONVOLUTIONS_1D_SSE_H
#define MATH_ALGORITHMS_NFFT_ATOMIC_CONVOLUTIONS_1D_SSE_H

#ifdef USE_SSE_ACCELERATION
#include <immintrin.h>
#endif

namespace MATH_ALGORITHMS
{
  namespace NFFT
  {
    template<int max_count, int count>
    struct nfft_atomic_convolution_sse
    {
      // maxcount is uneven !!
      inline static void execute(double* f, double* y, double* M){
        f[count] +=
          M[0+4*count]*y[0]
          +M[1+4*count]*y[1]
          +M[2+4*count]*y[2]
          +M[3+4*count]*y[3];

        //cout << M[0+4*count] << "\t" << M[1+4*count] << "\t" << M[2+4*count] << "\t" << M[3+4*count] << "\t" << endl;

        nfft_atomic_convolution_sse<max_count, count+1>::execute_sse(f, y, M);
      }

      // maxcount is uneven !!
      inline static void execute(double* f, double* y, double* A, double* B, double* C, double* D){
        f[count] +=
          A[count]*y[0]
          +B[count]*y[1]
          +C[count]*y[2]
          +D[count]*y[3];

        //cout << A[count] << "\t" << B[count] << "\t" << C[count] << "\t" << D[count] << "\t" << endl;

        nfft_atomic_convolution_sse<max_count, count+1>::execute_sse(f, y, A, B, C, D);
      }

#ifdef USE_SSE_ACCELERATION

      // y is vector of length 4
      // M is a  matrix (row-major), with 4 columns and max_count rows.
      inline static void execute_sse(double* f, double* y, double* M){

        __m128d F0, Y0, M0, M1;

        F0 = _mm_loadu_pd(f+count);

        { // col 0 and 1
          Y0 = _mm_load_pd(y+0);

          M0 = _mm_load_pd(M+4*count+0);
          M0 = _mm_mul_pd (M0,Y0);

          M1 = _mm_load_pd(M+4*count+4);
          M1 = _mm_mul_pd (M1,Y0);

          Y0 = _mm_hadd_pd(M0,M1);
          F0 = _mm_add_pd (F0, Y0);
        }

        { // col 2 and 3
          Y0 = _mm_load_pd(y+2);

          M0 = _mm_load_pd(M+4*count+2);
          M0 = _mm_mul_pd (M0,Y0);

          M1 = _mm_load_pd(M+4*count+6);
          M1 = _mm_mul_pd (M1,Y0);

          Y0 = _mm_hadd_pd(M0,M1);
          F0 = _mm_add_pd (F0, Y0);
        }

        _mm_storeu_pd(f+count, F0);

        nfft_atomic_convolution_sse<max_count, count+2>::execute_sse(f, y, M);
      }
#else
      inline static void execute_sse(double* f, double* y, double* M){
        f[count] +=
          M[0+4*count]*y[0]
          +M[1+4*count]*y[1]
          +M[2+4*count]*y[2]
          +M[3+4*count]*y[3];

        //cout << M[0+4*count] << "\t" << M[1+4*count] << "\t" << M[2+4*count] << "\t" << M[3+4*count] << "\t" << endl;

        nfft_atomic_convolution_sse<max_count, count+1>::execute_sse(f, y, M);
      }
#endif

      inline static void execute_sse(double* f, double* y, double* A, double* B, double* C, double* D){
        f[count] +=
          A[count]*y[0]
          +B[count]*y[1]
          +C[count]*y[2]
          +D[count]*y[3];

        //cout << A[count] << "\t" << B[count] << "\t" << C[count] << "\t" << D[count] << "\t" << endl;

        nfft_atomic_convolution_sse<max_count, count+1>::execute_sse(f, y, A, B, C, D);
      }

    };

    template<int max_count>
    struct nfft_atomic_convolution_sse<max_count, max_count>
    {
      template<typename scalar_type>
      inline static void execute_sse(scalar_type* f, scalar_type* y, scalar_type* M)
      {
        //cout << endl;
      }

      template<typename scalar_type>
      inline static void execute_sse(scalar_type* f, scalar_type* y, scalar_type* alpha, scalar_type* beta, scalar_type* gamma, scalar_type* delta)
      {
        //cout << endl;
      }
    };

  }

}

#endif
