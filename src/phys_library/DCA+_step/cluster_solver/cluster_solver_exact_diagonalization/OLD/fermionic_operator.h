//-*-C++-*-

#ifndef FERMIONIC_OPERATOR_H
#define FERMIONIC_OPERATOR_H

/*!
 *   \author Peter Staar
 */
template<typename b_dmn, typename s_dmn, typename r_dmn>
class fermionic_operator
{
  typedef fermionic_operator<b_dmn, s_dmn, r_dmn> this_type;
  
public:
  
  fermionic_operator();
  fermionic_operator(int b_ind, int s_ind, int r_ind);

  ~fermionic_operator();
  
  this_type& operator=(const this_type& rhs);

  //bool       operator<(const this_type& lhs, const this_type& rhs);

public:

  int linear_index;

  int sizes[3];
  int index[3];

  int Sz;
};

template<typename b_dmn, typename s_dmn, typename r_dmn>
fermionic_operator<b_dmn, s_dmn, r_dmn>::fermionic_operator()
{
  linear_index = -1;

  sizes[0] = b_dmn::dmn_size();
  sizes[1] = s_dmn::dmn_size();
  sizes[2] = r_dmn::dmn_size();

  for(int l=0; l<3; l++)
    index[l] = -1;
}

template<typename b_dmn, typename s_dmn, typename r_dmn>
fermionic_operator<b_dmn, s_dmn, r_dmn>::fermionic_operator(int b_ind, int s_ind, int r_ind)
{
  linear_index = (b_ind + sizes[0]*(s_ind + sizes[1]*r_ind)); 

  index[0] = b_ind;
  index[1] = s_ind;
  index[2] = r_ind;

  Sz = (s_ind==0)? -1 : 1;
}

template<typename b_dmn, typename s_dmn, typename r_dmn>
fermionic_operator<b_dmn, s_dmn, r_dmn>::~fermionic_operator()
{}

template<typename b_dmn, typename s_dmn, typename r_dmn>
fermionic_operator<b_dmn, s_dmn, r_dmn>& fermionic_operator<b_dmn, s_dmn, r_dmn>::operator=(fermionic_operator<b_dmn, s_dmn, r_dmn> const& rhs)
{
  linear_index = rhs.linear_index;

  index[0] = rhs.index[0];
  index[1] = rhs.index[1];
  index[2] = rhs.index[2];

  Sz = rhs.Sz;

  return (*this);
}

// template<typename b_dmn, typename s_dmn, typename r_dmn>
// bool fermionic_operator<b_dmn, s_dmn, r_dmn>::operator<(const fermionic_operator<b_dmn, s_dmn, r_dmn>& lhs,
// 							const fermionic_operator<b_dmn, s_dmn, r_dmn>& rhs)
// {
//   return (lhs.linear_index < rhs.linear_index);
// }

#endif
