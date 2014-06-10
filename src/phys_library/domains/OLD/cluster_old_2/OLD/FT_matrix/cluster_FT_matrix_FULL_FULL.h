//-*-C++-*-

/*
 *      Author: peterstaar
 */

#ifndef CLUSTER_FT_MATRIX_FULL_FULL_H_
#define CLUSTER_FT_MATRIX_FULL_FULL_H_

#include "cluster_FT_matrix.h"
#include <complex>
#include <math.h>

using namespace std;

template<class cluster_type>
class cluster_FT_matrix<cluster_type, FULL, FULL>
{
public:
  
  static std::complex<double>* get_matrix();
  static std::complex<double>* get_matrix_inverse();

private:

  static std::complex<double>* initialize();
  static std::complex<double>* initialize_inverse();

};

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, FULL, FULL>::get_matrix()
{
  static std::complex<double>* matrix = initialize();
  return matrix;
}

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, FULL, FULL>::get_matrix_inverse()
{
  static std::complex<double>* matrix = initialize_inverse();
  return matrix;
}

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, FULL, FULL>::initialize()
{
  // r --> k
  static std::complex<double>* matrix_ptr = new std::complex<double>[cluster_type::get_cluster_size()*cluster_type::get_cluster_size()];

  std::vector<double> r;
  std::vector<double> k;

  for(int i=0; i<cluster_type::get_cluster_size(); i++){
    for(int j=0; j<cluster_type::get_cluster_size(); j++){

      r = cluster_type::get_r_cluster()[i];
      k = cluster_type::get_k_cluster()[j];
      
      double rk = 0.;
      for(int l=0; l< int(r.size()); l++)
	  rk += r[l]*k[l];

      std::complex<double> c(cos(rk),sin(rk));
      matrix_ptr[j + cluster_type::get_cluster_size()*i] = c;
    }
  }

  return matrix_ptr;
}

template<class cluster_type>
std::complex<double>* cluster_FT_matrix<cluster_type, FULL, FULL>::initialize_inverse()
{
  // k --> r
  static std::complex<double>* matrix_ptr = new std::complex<double>[cluster_type::get_cluster_size()*cluster_type::get_cluster_size()];

  std::vector<double> r;
  std::vector<double> k;

  for(int i=0; i<cluster_type::get_cluster_size(); i++){
    for(int j=0; j<cluster_type::get_cluster_size(); j++){

      r = cluster_type::get_r_cluster()[i];
      k = cluster_type::get_k_cluster()[j];
      
      double rk = 0.;
      for(int l=0; l<int(r.size()); l++)
	  rk += r[l]*k[l];

      std::complex<double> c_inv(cos(rk),-sin(rk));
      matrix_ptr[i + cluster_type::get_cluster_size()*j] = c_inv/double(cluster_type::get_cluster_size());
    }
  }

  return matrix_ptr;
}

#endif
