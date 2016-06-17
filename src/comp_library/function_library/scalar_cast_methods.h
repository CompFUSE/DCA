//-*-C++-*-
#ifndef SCALAR_CAST_METHODS_H
#define SCALAR_CAST_METHODS_H
#include <complex>
/*!
 *  \author: peterstaar
 *    
 *  cast to a real floating point number
 */
template<typename scalartype>
struct do_cast
{
  template<typename scalartype_2>
  static scalartype execute(scalartype_2 value);

  template<typename scalartype_2>
  static scalartype execute(std::complex<scalartype_2> value);
};


template<typename scalartype>
template<typename scalartype_2>
scalartype do_cast<scalartype>::execute(scalartype_2 value)
{
  return static_cast<scalartype>(value);
}

template<typename scalartype>
template<typename scalartype_2>
scalartype do_cast<scalartype>::execute(std::complex<scalartype_2> value)
{
  return static_cast<scalartype>(real(value) );
}

/*!
 *  \author: Peter Staar
 *
 *  cast to a complex floating point number
 */
template<typename scalartype>
struct do_cast<std::complex<scalartype> >
{
  template<typename scalartype_2>
  static inline std::complex<scalartype> execute(scalartype_2 value);

  template<typename scalartype_2>
  static inline std::complex<scalartype> execute(std::complex<scalartype_2> value);
};

template<typename scalartype>
template<typename scalartype_2>
std::complex<scalartype> do_cast<std::complex<scalartype> >::execute(scalartype_2 value)
{
  return static_cast<std::complex<scalartype> >(value);
}

template<typename scalartype>
template<typename scalartype_2>
std::complex<scalartype> do_cast<std::complex<scalartype> >::execute(std::complex<scalartype_2> value)
{
  return static_cast<std::complex<scalartype> >(value);
}


#endif
