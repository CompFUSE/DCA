//-*-C++-*-

/*
 *      Author: peterstaar
 */

#ifndef CLUSTER_FT_MATRIX_IRREDUCIBLE_FULL_H_
#define CLUSTER_FT_MATRIX_IRREDUCIBLE_FULL_H_

#include "cluster_FT_matrix.h"
#include <complex>
#include <math.h>

using namespace std;

template<class cluster_type>
class cluster_FT_matrix<cluster_type, IRREDUCIBLE, FULL>
{
public:
  
  static std::complex<double>* get_matrix();
  static std::complex<double>* get_matrix_inverse();

private:

  static std::complex<double>*  initialize();
  static std::complex<double>*  initialize_inverse();

};

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, IRREDUCIBLE, FULL>::get_matrix()
{
  static std::complex<double>* matrix = initialize();
  return matrix;
}

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, IRREDUCIBLE, FULL>::get_matrix_inverse()
{
  static std::complex<double>* matrix = initialize_inverse();
  return matrix;
}

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, IRREDUCIBLE, FULL>::initialize()
{
  assert(false);

  static std::complex<double>* matrix_ptr = new std::complex<double>[cluster_type::get_cluster_size()*cluster_type::get_irreducible_cluster_size()];
  memset(matrix_ptr,0.,sizeof(std::complex<double>)*cluster_type::get_irreducible_cluster_size()*cluster_type::get_cluster_size());

  std::vector<double> r;
  std::vector<double> k;

  for(int i=0; i<cluster_type::get_cluster_size(); i++){
    for(int j=0; j<cluster_type::get_cluster_size(); j++){

      r = cluster_type::get_r_cluster()[i];
      k = cluster_type::get_k_cluster()[j];
      
      double rk = 0.;
      for(int l=0; l<int(r.size()); l++)
	  rk += r[l]*k[l];

      int i_ind = cluster_type::get_conversion_table()[i];
      int matrix_ind_i = cluster_type::get_irreducible_index_from_full_index(i);

      std::complex<double> c(cos(rk),sin(rk));
      matrix_ptr[j + cluster_type::get_cluster_size()*matrix_ind_i] += c/double(cluster_type::get_weights()[i_ind]);
    }
  }

  return matrix_ptr;
}

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, IRREDUCIBLE, FULL>::initialize_inverse()
{
  assert(false);

  // k(IRREDUCIBLE) --> r(FULL)
  static std::complex<double>* matrix_ptr = new std::complex<double>[cluster_type::get_cluster_size()*cluster_type::get_irreducible_cluster_size()];
  memset(matrix_ptr,0.,sizeof(std::complex<double>)*cluster_type::get_irreducible_cluster_size()*cluster_type::get_cluster_size());

  std::vector<double> r;
  std::vector<double> k;

  for(int i=0; i<cluster_type::get_cluster_size(); i++){
    for(int j=0; j<cluster_type::get_cluster_size(); j++){

      r = cluster_type::get_r_cluster()[i];
      k = cluster_type::get_k_cluster()[j];
      
      double rk = 0.;
      for(int l=0; l<int(r.size()); l++)
	  rk += r[l]*k[l];

      //int j_ind = cluster_type::get_conversion_table()[j];
      int matrix_ind_j = cluster_type::get_irreducible_index_from_full_index(j);

      std::complex<double> c_inv(cos(rk),-sin(rk));
      matrix_ptr[i + cluster_type::get_cluster_size()*matrix_ind_j] += c_inv/double(cluster_type::get_cluster_size());
    }
  }

  return matrix_ptr;
}

#endif
