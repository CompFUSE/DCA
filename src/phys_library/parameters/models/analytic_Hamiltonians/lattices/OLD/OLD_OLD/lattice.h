//-*-C++-*-

/*
 *      Author: peter staar
 */

/*
#ifndef LATTICE_H_
#define LATTICE_H_

// lattice-types
#include "2D_square_A.h"
#include "2D_triangular_A.h"
#include "2D_parallellogrammatic_A.h"

// each lattice-type will define the initialization-functions
template<typename lattice_type>
class lattice
{
public:

  const static int DIMENSION = lattice_type::DIMENSION;
  static int BANDS;

  //  const static double Fermi_energy = 0.;
  //  const static double density      = 1.;

  static double* get_r_DCA_basis();
  static double* get_k_DCA_basis();

  static double* get_r_LDA_basis();
  static double* get_k_LDA_basis();

  static std::complex<double> get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2);

  static void                 symmetrize_Hamiltonian(std::complex<double>* H_matrix);

  static std::vector<int>     get_interacting_bands();

};

template<typename lattice_type>
int lattice<lattice_type>::BANDS = lattice_type::BANDS;

template<typename lattice_type>
double* lattice<lattice_type>::get_r_DCA_basis()
{
  static double* r_DCA = lattice_type::get_r_DCA_basis();
  return r_DCA;
}

template<typename lattice_type>
double* lattice<lattice_type>::get_k_DCA_basis()
{
  static double* k_DCA = lattice_type::get_k_DCA_basis(); 
  return k_DCA;
}

template<typename lattice_type>
double* lattice<lattice_type>::get_r_LDA_basis()
{
  static double* r_LDA = lattice_type::get_r_LDA_basis();
  return r_LDA;
}

template<typename lattice_type>
double* lattice<lattice_type>::get_k_LDA_basis()
{
  static double* k_LDA = lattice_type::get_k_LDA_basis();
  return k_LDA;
}

template<typename lattice_type>
std::complex<double> lattice<lattice_type>::get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2)
{
  return lattice_type::get_LDA_Hamiltonians(t, k, b1, s1, b2, s2);
}

template<typename lattice_type>
void lattice<lattice_type>::symmetrize_Hamiltonian(std::complex<double>* H_matrix)
{
  lattice_type::symmetrize_Hamiltonian(H_matrix);
}




#endif
*/
