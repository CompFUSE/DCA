//-*-C++-*-

#ifndef SCALAR_CAST_METHODS_H
#define SCALAR_CAST_METHODS_H

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

/*
 *    cast to an integer number
 */
// template<>
// struct do_cast<int>
// {

//   //   template<typename scalartype_1>
//   //   static inline int execute(scalartype_1 value_1);
// };

/*
  template<typename scalartype_1>
  int do_cast<int>::execute(scalartype_1 value_1)
  {
  throw std::logic_error(__FUNCTION__);
  return -1;
  }


  template<>
  int do_cast<int>::execute(e_spin_states_type e_spin)
  {
  switch(e_spin)
  {
  case e_DN:
  return 0;
  break;

  case e_UP:
  return 1;
  break;

  default:
  throw std::logic_error(__FUNCTION__);
  }
  }


  template<>
  int do_cast<int>::execute(HS_spin_states_type HS_spin)
  {
  switch(HS_spin)
  {
  case HS_DN:
  return 0;
  break;

  case HS_ZERO:
  return 1;
  break;

  case HS_UP:
  return 2;
  break;

  default:
  throw std::logic_error(__FUNCTION__);
  }
  }

  template<>
  int do_cast<int>::execute(HS_field_sign_type HS_field_sign)
  {
  switch(HS_field_sign)
  {
  case HS_FIELD_DN:
  return 0;
  break;

  case HS_FIELD_UP:
  return 1;
  break;

  default:
  throw std::logic_error(__FUNCTION__);
  }
  }
*/

#endif
