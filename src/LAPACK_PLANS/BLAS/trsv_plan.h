//-*-C++-*-

/*
 * triangular_solve_plan.h
 *
 *      Author: peter staar
 */

#ifndef TRIANGULAR_SOLVE_PLAN_H_
#define TRIANGULAR_SOLVE_PLAN_H_

/*
template<typename scalartype>
class triangular_solve_plan
{
public:

  static void execute(char uplo, char trans, char diag, int n, scalartype* A, int LDA, scalartype* X, int incx);

};

template<typename scalartype>
void triangular_solve_plan<scalartype>::execute(char uplo, char trans, char diag, int n, scalartype* A, int LDA, scalartype* X, int incx)
{
  throw std::logic_error(__FUNCTION__);
}

template<>
void triangular_solve_plan<float>::execute(char uplo, char trans, char diag, int n, float* A, int LDA, float* X, int incx)
{
  BLAS::strsv_(&uplo, &trans, &diag, &n, A, &LDA , X, &incx);
}

template<>
void triangular_solve_plan<double>::execute(char uplo, char trans, char diag, int n, double* A, int LDA, double* X, int incx)
{
  BLAS::dtrsv_(&uplo, &trans, &diag, &n, A, &LDA , X, &incx);
}

template<>
void triangular_solve_plan<std::complex<float> >::execute(char uplo, char trans, char diag, int n, std::complex<float>* A, int LDA, std::complex<float>* X, int incx)
{
  BLAS::ctrsv_(&uplo, &trans, &diag, &n, A, &LDA , X, &incx);
}

template<>
void triangular_solve_plan<std::complex<double> >::execute(char uplo, char trans, char diag, int n, std::complex<double>* A, int LDA, std::complex<double>* X, int incx)
{
  BLAS::ztrsv_(&uplo, &trans, &diag, &n, A, &LDA , X, &incx);
}
*/

#endif
