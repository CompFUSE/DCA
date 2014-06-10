//-*-C++-*-

#ifndef CLUSTER_HARMONIC_R_BASIS_FUNCTIONS_H
#define CLUSTER_HARMONIC_R_BASIS_FUNCTIONS_H

/*!
 *  \ingroup CLUSTER
 *
 *  \author  Peter Staar
 */
template<class basis_cluster_type>
class cluster_harmonics_r_basis_functions
{
  const static int DIMENSION = basis_cluster_type::DIMENSION;

  typedef r_cluster<FULL,basis_cluster_type>   source_r_cluster_cluster_type;
  typedef dmn_0<source_r_cluster_cluster_type> source_r_dmn_t;

  typedef k_cluster<FULL,basis_cluster_type>   source_k_cluster_cluster_type;
  typedef dmn_0<source_k_cluster_cluster_type> source_k_dmn_t;

public:

  typedef std::vector<std::vector<double> >                       element_type;
  typedef cluster_harmonics_r_basis_functions<basis_cluster_type> this_type;
  
public:

  static int                        get_size();
  static std::vector<element_type>& get_elements();

private:

  static std::vector<element_type> initialize();
};

template<class basis_cluster_type>
int cluster_harmonics_r_basis_functions<basis_cluster_type>::get_size()
{
  return source_r_dmn_t::dmn_size();
}

template<class basis_cluster_type>
std::vector<typename cluster_harmonics_r_basis_functions<basis_cluster_type>::element_type>&
cluster_harmonics_r_basis_functions<basis_cluster_type>::get_elements()
{
  static std::vector<element_type> elements = initialize();
  return elements;
}

template<class basis_cluster_type>
std::vector<typename cluster_harmonics_r_basis_functions<basis_cluster_type>::element_type>
cluster_harmonics_r_basis_functions<basis_cluster_type>::initialize()
{
  std::vector<element_type> Phi_r(source_r_dmn_t::dmn_size());

  for(int phi_ind=0; phi_ind<source_r_dmn_t::dmn_size(); phi_ind++){
    Phi_r[phi_ind] = cluster_harmonics<basis_cluster_type>::get_basis_function_in_real_space(phi_ind);
  }

  return Phi_r;
}





#endif
