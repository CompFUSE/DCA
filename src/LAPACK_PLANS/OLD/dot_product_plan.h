//-*-C++-*-

/*
 * dot_product_plan.h
 *
 *      Author: peter staar
 */




#ifndef DOT_PRODUCT_PLAN_H_
#define DOT_PRODUCT_PLAN_H_


template<typename scalartype>
class dot_product_plan
{

 public:

  dot_product_plan(int n);
  ~dot_product_plan();

  scalartype execute_plan();
  static scalartype execute_plan(int N, scalartype* X, int INCX, scalartype* Y, int INCY);

  int N;
  int INCX;
  int INCY;
  
  scalartype* X;
  scalartype* Y;
};

extern "C"              float   sdot_(const int *N, const float*                X, const int* INCX, const float*                Y,  const int* INCY);
extern "C"              double  ddot_(const int *N, const double*               X, const int* INCX, const double*               Y,  const int* INCY);
extern "C" std::complex<float > cdotc_(const int *N, const std::complex<float>*  X, const int* INCX, const std::complex<float>*  Y,  const int* INCY);
extern "C" std::complex<double> zdotc_(const int *N, const std::complex<double>* X, const int* INCX, const std::complex<double>* Y,  const int* INCY);

template<typename scalartype>
dot_product_plan<scalartype>::dot_product_plan(int n):
  N(n),
  INCX(1),
  INCY(1)
{}

template<typename scalartype>
dot_product_plan<scalartype>::~dot_product_plan()
{}

template<typename scalartype>
scalartype dot_product_plan<scalartype>::execute_plan()
{
  throw std::logic_error(__FUNCTION__);
  return 0;
}

template<>
float dot_product_plan<float>::execute_plan()
{
  return sdot_(&N, X, &INCX, Y, &INCY);
}

template<>
double dot_product_plan<double>::execute_plan()
{
  return ddot_(&N, X, &INCX, Y, &INCY);
}

template<>
std::complex<float> dot_product_plan<std::complex<float> >::execute_plan()
{
  return cdotc_(&N, X, &INCX, Y, &INCY);
}

template<>
std::complex<double> dot_product_plan<std::complex<double> >::execute_plan()
{
  return zdotc_(&N, X, &INCX, Y, &INCY);
}

template<typename scalartype>
scalartype dot_product_plan<scalartype>::execute_plan(int N, scalartype* X, int INCX, scalartype* Y, int INCY)
{
  throw std::logic_error(__FUNCTION__);
  return 0;
}

template<>
float dot_product_plan<float>::execute_plan(int N, float* X, int INCX, float* Y, int INCY)
{
  return sdot_(&N, X, &INCX, Y, &INCY);
}

template<>
double dot_product_plan<double>::execute_plan(int N, double* X, int INCX, double* Y, int INCY)
{
  return ddot_(&N, X, &INCX, Y, &INCY);
}

template<>
std::complex<float> dot_product_plan<std::complex<float> >::execute_plan(int N, std::complex<float>* X, int INCX, std::complex<float>* Y, int INCY)
{
  return cdotc_(&N, X, &INCX, Y, &INCY);
}

template<>
std::complex<double> dot_product_plan<std::complex<double> >::execute_plan(int N, std::complex<double>* X, int INCX, std::complex<double>* Y, int INCY)
{
  return zdotc_(&N, X, &INCX, Y, &INCY);
}

#endif
