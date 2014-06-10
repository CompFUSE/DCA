//-*-C++-*-

#ifndef LEGENDRE_DOMAIN_H
#define LEGENDRE_DOMAIN_H

//#include "legendre_basis_functions.h"

/*!
 *      Author: bart ydens
 */
/*
template<typename time_domain_type, legendre_representation_type L_rep_type=SINGLE_PARTICLE_QUANTITY>
class legendre_domain 
{
public:

  typedef int                                        element_type;
  typedef legendre_basis_functions<time_domain_type> basis_function_t;

public:

  static int&              get_size();
  static std::vector<int>& get_elements();

  template<typename parameters_type>
  static void initialize(parameters_type& parameters);

  static std::vector<int>& initialize_elements();

  template<class stream_type>
  static void to_JSON(stream_type& ss);
};

template<typename time_domain_type, legendre_representation_type L_rep_type>
int& legendre_domain<time_domain_type, L_rep_type>::get_size()  
{
  static int size =-1;
  return size;
}

template<typename time_domain_type, legendre_representation_type L_rep_type>
std::vector<int>& legendre_domain<time_domain_type, L_rep_type>::get_elements()
{
  static std::vector<int>& v = initialize_elements();
  return v;
}

template<typename time_domain_type, legendre_representation_type L_rep_type>
template<typename parameters_type>
void legendre_domain<time_domain_type, L_rep_type>::initialize(parameters_type& parameters)
{
  get_size() = parameters.get_nb_of_legendre_coefficients_single_particle();
}

template<typename time_domain_type, legendre_representation_type L_rep_type>
std::vector<int>& legendre_domain<time_domain_type, L_rep_type>::initialize_elements()
{
  static std::vector<int> v(get_size());

  for(int i=0; i<get_size(); i++)
    v[i] = i;

  return v;
}

template<typename time_domain_type, legendre_representation_type L_rep_type>
template<class stream_type>
void legendre_domain<time_domain_type, L_rep_type>::to_JSON(stream_type& ss)
{
  ss << "\"legendre_domain\" : [\n";
    
  for(int i=0; i<get_size(); i++)
    if(i == get_size()-1)
      ss << get_elements()[i] << "\n";
    else
      ss << get_elements()[i] << ",\n";
  
  ss << "]\n";
}
*/

#endif
