//-*-C++-*-

/*
 *      Author: peterstaar
 */

#ifndef CLUSTER_FT_MATRIX_IRREDUCIBLE_IRREDUCIBLE_H_
#define CLUSTER_FT_MATRIX_IRREDUCIBLE_IRREDUCIBLE_H_

#include "cluster_FT_matrix.h"
#include <complex>
#include <math.h>

using namespace std;

template<class cluster_type>
class cluster_FT_matrix<cluster_type, IRREDUCIBLE, IRREDUCIBLE>
{
public:
  
  static std::complex<double>* get_matrix();
  static std::complex<double>* get_matrix_inverse();

private:

  static std::complex<double>*  initialize();
  static std::complex<double>* initialize_inverse();

};

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, IRREDUCIBLE, IRREDUCIBLE>::get_matrix()
{
  static std::complex<double>* matrix = initialize();
  return matrix;
}

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, IRREDUCIBLE, IRREDUCIBLE>::get_matrix_inverse()
{
  static std::complex<double>* matrix = initialize_inverse();
  return matrix;
}

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, IRREDUCIBLE, IRREDUCIBLE>::initialize()
{
  // r --> k
  static std::complex<double>* matrix_ptr = new std::complex<double>[cluster_type::get_irreducible_cluster_size()*cluster_type::get_irreducible_cluster_size()];
  memset(matrix_ptr,0.,sizeof(std::complex<double>)*cluster_type::get_irreducible_cluster_size()*cluster_type::get_irreducible_cluster_size());

  std::vector<double> r;
  std::vector<double> k;

  for(int i=0; i<cluster_type::get_cluster_size(); i++){
    for(int j=0; j<cluster_type::get_cluster_size(); j++){

      r = cluster_type::get_r_cluster()[i];
      k = cluster_type::get_k_cluster()[j];
      
      double rk = 0.;
      for(int l=0; l<int(r.size()); l++)
	  rk += r[l]*k[l];

      //int i_ind = cluster_type::get_conversion_table()[i];
      int j_ind = cluster_type::get_conversion_table()[j];

      int matrix_ind_i = cluster_type::get_irreducible_index_from_full_index(i);
      int matrix_ind_j = cluster_type::get_irreducible_index_from_full_index(j);
      
      std::complex<double> c(cos(rk),sin(rk));
      matrix_ptr[matrix_ind_j + cluster_type::get_irreducible_cluster_size()*matrix_ind_i] += c/double(cluster_type::get_weights()[j_ind]);
    }
  }

  return matrix_ptr;
}

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, IRREDUCIBLE, IRREDUCIBLE>::initialize_inverse()
{
  // k --> r
  static std::complex<double>* matrix_ptr = new std::complex<double>[cluster_type::get_irreducible_cluster_size()*cluster_type::get_irreducible_cluster_size()];
  memset(matrix_ptr,0.,sizeof(std::complex<double>)*cluster_type::get_irreducible_cluster_size()*cluster_type::get_irreducible_cluster_size());

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
      //int j_ind = cluster_type::get_conversion_table()[j];

      int matrix_ind_i = cluster_type::get_irreducible_index_from_full_index(i);
      int matrix_ind_j = cluster_type::get_irreducible_index_from_full_index(j);

      std::complex<double> c_inv(cos(rk),-sin(rk));
      matrix_ptr[matrix_ind_i + cluster_type::get_irreducible_cluster_size()*matrix_ind_j] += c_inv/double(cluster_type::get_cluster_size())/double(cluster_type::get_weights()[i_ind]);
    }
  }

  return matrix_ptr;
}

#endif
