//-*-C++-*-

#ifndef BLAS_GEMM_H
#define BLAS_GEMM_H

/*!
 * matrix_matrix_multiplication_plan.h
 *
 *      Author: peter staar
 */
template<typename scalartype>
void blas_gemm_t(char& TRANSA, char& TRANSB, int& M, int& N, int& K, scalartype& alpha,
		 scalartype* A, int& LDA, scalartype* B, int& LDB, scalartype& beta,
		 scalartype* C, int& LDC)
{
  cout << __PRETTY_FUNCTION__ << endl;
  throw std::logic_error(__FUNCTION__);
}

template<>
void blas_gemm_t<float>(char& TRANSA, char& TRANSB, int& M, int& N, int& K, float& alpha,
			float* A, int& LDA, float* B, int& LDB, float& beta,
			float* C, int& LDC)
{
  BLAS::sgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}
  
template<>
void blas_gemm_t<double>(char& TRANSA, char& TRANSB, int& M, int& N, int& K, double& alpha,
			 double* A, int& LDA, double* B, int& LDB, double& beta,
			 double* C, int& LDC)
{
  BLAS::dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}
  
template<>
void blas_gemm_t<std::complex<float> >(char& TRANSA, char& TRANSB, int& M, int& N, int& K, std::complex<float>& alpha,
				       std::complex<float>* A, int& LDA, std::complex<float>* B, int& LDB, std::complex<float>& beta,
				       std::complex<float>* C, int& LDC)
{
  BLAS::cgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}
  
template<>
void blas_gemm_t<std::complex<double> >(char& TRANSA, char& TRANSB, int& M, int& N, int& K, std::complex<double>& alpha,
					std::complex<double>* A, int& LDA, std::complex<double>* B, int& LDB, std::complex<double>& beta,
					std::complex<double>* C, int& LDC)
{
  BLAS::zgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);
}

#endif
