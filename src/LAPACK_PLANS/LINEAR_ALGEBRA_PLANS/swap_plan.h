//-*-C++-*-

#ifndef SWAP_PLAN_H
#define SWAP_PLAN_H

/*!
 *  \author Peter Staar
 */
class swap_plan
{
public:

  template<typename scalar_type>
  static void execute(int length, scalar_type* a, int inc_a, scalar_type* b, int inc_b);
};

/*
template<typename scalar_type>
inline void swap_plan::execute(int length, scalar_type* a, int inc_a, scalar_type* b, int inc_b)
{
  throw std::logic_error(__FUNCTION__);
}

template<>
inline void swap_plan::execute(int length, float* a, int inc_a, float* b, int inc_b)
{
  BLAS::sswap_(&length, a, &inc_a, b, &inc_b);
}

template<>
inline void swap_plan::execute(int length, double* a, int inc_a, double* b, int inc_b)
{
  BLAS::dswap_(&length, a, &inc_a, b, &inc_b);
}

template<>
inline void swap_plan::execute(int length, std::complex<float>* a, int inc_a, std::complex<float>* b, int inc_b)
{
  BLAS::cswap_(&length, a, &inc_a, b, &inc_b);
}

template<>
inline void swap_plan::execute(int length, std::complex<double>* a, int inc_a, std::complex<double>* b, int inc_b)
{
  BLAS::zswap_(&length, a, &inc_a, b, &inc_b);
}
*/


#endif
