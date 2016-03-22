//-*-C++-*-

#ifndef SCALE_PLAN_H
#define SCALE_PLAN_H

/*!
 *  \author Peter Staar
 *
 *  \brief scales vector 'a' with factor f.
 */
/*
class scale_plan
{
public:

  template<typename scalar_type>
    static void execute(int length, scalar_type f, scalar_type* a, int inc_a);
};

template<typename scalar_type>
inline void scale_plan::execute(int length, scalar_type f, scalar_type* a, int inc_a)
{
  throw std::logic_error(__FUNCTION__);
}

template<>
inline void scale_plan::execute(int length, float f, float* a, int inc_a)
{
  BLAS::sscal_(&length, &f, a, &inc_a);
}

template<>
inline void scale_plan::execute(int length, double f, double* a, int inc_a)
{
  BLAS::dscal_(&length, &f, a, &inc_a);
}

template<>
inline void scale_plan::execute(int length, std::complex<float> f, std::complex<float>* a, int inc_a)
{
  BLAS::cscal_(&length, &f, a, &inc_a);
}

template<>
inline void scale_plan::execute(int length, std::complex<double> f, std::complex<double>* a, int inc_a)
{
  BLAS::zscal_(&length, &f, a, &inc_a);
}
*/

#endif
