//-*-C++-*-

#ifndef CLUSTER_HARMONICS_H
#define CLUSTER_HARMONICS_H

/*!
 *  \ingroup CLUSTER
 *
 *  \author  Peter Staar
 */
template<class basis_cluster_type>
class cluster_harmonics
{
  typedef r_cluster<FULL,basis_cluster_type>   source_r_cluster_cluster_type;
  typedef dmn_0<source_r_cluster_cluster_type> source_r_dmn_t;

  const static int DIMENSION = basis_cluster_type::DIMENSION;

public:

  typedef std::vector<double>                   element_type;
  typedef cluster_harmonics<basis_cluster_type> this_type;
  
public:

  static int                  get_size();
  static std::vector<double>& get_elements();

  static std::vector<std::vector<double> > get_basis_function_in_real_space(int n);
  static std::complex<double>              evaluate_at_k(int n, std::vector<double> k);

private:

  static FUNC_LIB::function<std::vector<std::vector<double> >, source_r_dmn_t >& initialize();
};

template<class basis_cluster_type>
int cluster_harmonics<basis_cluster_type>::get_size()
{
  return basis_cluster_type::get_cluster_size();
}

template<class basis_cluster_type>
std::vector<double>& cluster_harmonics<basis_cluster_type>::get_elements()
{
  return source_r_cluster_cluster_type::get_elements();
}

template<class basis_cluster_type>
std::complex<double> cluster_harmonics<basis_cluster_type>::evaluate_at_k(int n, std::vector<double> k)
{
  assert(n>=0 && n<get_size());

  std::complex<double>              phase_factor    = 0;
  std::vector<std::vector<double> > delta_functions = get_basis_function_in_real_space(n);

  for(size_t r_ind=0; r_ind<delta_functions.size(); r_ind++){
	      
    std::vector<double> R_vec = delta_functions[r_ind];

    double dt_prdct = VECTOR_OPERATIONS::DOT_PRODUCT(R_vec, k);
    
    phase_factor += std::complex<double>(cos(dt_prdct), sin(dt_prdct));
  }

//   if(n==0){
//     cout << phase_factor << endl;
    
//     if(abs(phase_factor-1.)>1.e-6)
//       throw std::logic_error(__FUNCTION__);
//   }

  return phase_factor/double(delta_functions.size());
}

template<class basis_cluster_type>
std::vector<std::vector<double> > cluster_harmonics<basis_cluster_type>::get_basis_function_in_real_space(int n)
{
  static FUNC_LIB::function<std::vector<std::vector<double> >, source_r_dmn_t >& centered_r_cluster = initialize();
  return centered_r_cluster(n);
}

template<class basis_cluster_type>
FUNC_LIB::function<std::vector<std::vector<double> >, typename cluster_harmonics<basis_cluster_type>::source_r_dmn_t >&
cluster_harmonics<basis_cluster_type>::initialize()
{
  static FUNC_LIB::function<std::vector<std::vector<double> >, source_r_dmn_t > centered_r_cluster("centered_r_cluster");

  for(int R_ind=0; R_ind<source_r_dmn_t::dmn_size(); R_ind++){
    
    centered_r_cluster(R_ind) = source_r_cluster_cluster_type::find_equivalent_vectors_with_minimal_distance_to_origin(source_r_dmn_t::get_elements()[R_ind]).second;
    
    int size = centered_r_cluster(R_ind).size();

    for(int r_ind=0; r_ind<size; r_ind++){
      std::vector<double> min_r = centered_r_cluster(R_ind)[r_ind];
      
      for(size_t l=0; l<min_r.size(); l++)
	min_r[l] *= -1.;
      
      centered_r_cluster(R_ind).push_back(min_r);
    }

    sort(centered_r_cluster(R_ind).begin(), centered_r_cluster(R_ind).end(), &VECTOR_OPERATIONS::IS_LARGER_VECTOR);
    int vec_size = unique(centered_r_cluster(R_ind).begin(), centered_r_cluster(R_ind).end(), &SAME_VECTOR) - centered_r_cluster(R_ind).begin();
    
    centered_r_cluster(R_ind).resize(vec_size);
  }

  return centered_r_cluster;
}

template<class k_cluster_type>
class cluster_harmonics<dmn_0<k_cluster_type> >
{
  typedef typename k_cluster_type::base_cluster base_cluster_type;

public:

  static std::complex<double> evaluate(int n, std::vector<double> k)
  {
    return cluster_harmonics<base_cluster_type>::evaluate_at_k(n, k);
  }
};

#endif
