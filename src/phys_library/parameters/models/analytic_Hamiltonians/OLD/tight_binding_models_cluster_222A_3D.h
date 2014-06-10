//-*-C++-*-

/*
 *      Author: peter staar
 */

/*
#ifndef tight_binding_models_cluster_222A_3D_H_
#define tight_binding_models_cluster_222A_3D_H_

//#include "tight_binding_models.h"

//#include "No_symmetry.h"
//#include "Oh.h"

template<>
class tight_binding_models<cluster_222A_3D> {

 public:

  typedef Oh                DCA_point_group;
  typedef Null_symmetry_3D  LDA_point_group;

  static double Fermi_energy;

  const static int DIMENSION = 3;
  const static int BANDS     = 1;

  static int get_DCA_size();
  static int get_LDA_size();

  static std::vector<int>& DCA_grid_size();
  static std::vector<int>& LDA_grid_size();

  static double* get_r_DCA_basis();
  static double* get_k_DCA_basis();

  static double* get_r_LDA_basis();
  static double* get_k_LDA_basis();

  static std::complex<double> get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2);

  static void symmetrize_Hamiltonian(std::complex<double>* H_matrix);

  static double interaction_Hamiltonian(int nu_i, int nu_j, int delta_r=0);
  static double alpha_field            (int nu_i, int nu_j, int delta_r=0);

  template<class parameters_type>
  static void   initialize_interaction_Hamiltonian(parameters_type& parameters);

  static std::vector<int> get_interacting_bands();

private:
  static double U;
};

double tight_binding_models<cluster_222A_3D>::U = 0;

double tight_binding_models<cluster_222A_3D>::Fermi_energy = 0.;

int tight_binding_models<cluster_222A_3D>::get_DCA_size()
{
  return DCA_grid_size()[0]*DCA_grid_size()[1]*DCA_grid_size()[2];
}

int tight_binding_models<cluster_222A_3D>::get_LDA_size()
{
  return LDA_grid_size()[0]*LDA_grid_size()[1]*LDA_grid_size()[2];
}

int get_LDA_size();

std::vector<int>& tight_binding_models<cluster_222A_3D>::DCA_grid_size()
{
  static std::vector<int> v(3);
  v[0] = 2;
  v[1] = 2;
  v[2] = 2;

  return v;
}

std::vector<int>& tight_binding_models<cluster_222A_3D>::LDA_grid_size()
{
  static std::vector<int> v(3);
  v[0] = 20;
  v[1] = 20;
  v[2] = 20;

  return v;
}

double* tight_binding_models<cluster_222A_3D>::get_r_DCA_basis()
{
  static double* r_DCA = new double[9];
  
  r_DCA[0] = 1.;  r_DCA[1] = 0.; r_DCA[2] = 0.;
  r_DCA[3] = 0.;  r_DCA[4] = 1.; r_DCA[5] = 0.;
  r_DCA[6] = 0.;  r_DCA[7] = 0.; r_DCA[8] = 1.;

  return r_DCA;
}

double* tight_binding_models<cluster_222A_3D>::get_k_DCA_basis()
{
   static double* k_DCA = new double[9];
  
  k_DCA[0] = 2.*M_PI;  k_DCA[1] = 0.;      k_DCA[2] = 0.;
  k_DCA[3] = 0.;       k_DCA[4] = 2.*M_PI; k_DCA[5] = 0.;
  k_DCA[6] = 0.;       k_DCA[7] = 0.;      k_DCA[8] = 2.*M_PI;

  return k_DCA;
}

double* tight_binding_models<cluster_222A_3D>::get_r_LDA_basis()
{
  static double* r_LDA = new double[9];
  
  r_LDA[0] = 1.;  r_LDA[1] = 0.; r_LDA[2] = 0.;
  r_LDA[3] = 0.;  r_LDA[4] = 1.; r_LDA[5] = 0.;
  r_LDA[6] = 0.;  r_LDA[7] = 0.; r_LDA[8] = 1.;

  return r_LDA;
}

double* tight_binding_models<cluster_222A_3D>::get_k_LDA_basis()
{
   static double* k_LDA = new double[9];
  
  k_LDA[0] = 2.*M_PI;  k_LDA[1] = 0.;      k_LDA[2] = 0.;
  k_LDA[3] = 0.;       k_LDA[4] = 2.*M_PI; k_LDA[5] = 0.;
  k_LDA[6] = 0.;       k_LDA[7] = 0.;      k_LDA[8] = 2.*M_PI;

  return k_LDA;
}

std::complex<double> tight_binding_models<cluster_222A_3D>::get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2)
{
  std::complex<double> H_LDA = 0.;

  if( (b1 == 0 && b2 == 0) && ((s1 == 0 && s2 == 0) || (s1 == 1 && s2 == 1)))
    H_LDA = -2.* t * ( cos(k[0]) + cos(k[1]) + cos(k[2]) );

  return H_LDA;
}

void tight_binding_models<cluster_222A_3D>::symmetrize_Hamiltonian(std::complex<double>* H_matrix)
{
  //const static int matrix_size = (2*BANDS)*(2*BANDS);
  const static int matrix_dim  = (2*BANDS);

  // symmetrize the e_spin up and dn part
  std::complex<double> tmp = (H_matrix[0+matrix_dim*0] + H_matrix[1+matrix_dim*1])/double(2); 

  H_matrix[0+matrix_dim*0] = tmp;
  H_matrix[1+matrix_dim*0] = 0;
  H_matrix[0+matrix_dim*1] = 0;
  H_matrix[1+matrix_dim*1] = tmp;
}

double tight_binding_models<cluster_222A_3D>::interaction_Hamiltonian(int nu_i, int nu_j, int delta_r)
{
  if( (nu_i==nu_j) && (delta_r==0) )
    return 0.;
  else
    return U;
}

template<class parameters_type>
void tight_binding_models<cluster_222A_3D>::initialize_interaction_Hamiltonian(parameters_type& parameters)
{
  U = parameters.get_U_hubbard();
}


std::vector<int> tight_binding_models<cluster_222A_3D>::get_interacting_bands()
{
  std::vector<int> v;
  v.push_back(0);
  return v;
}

#endif
*/
