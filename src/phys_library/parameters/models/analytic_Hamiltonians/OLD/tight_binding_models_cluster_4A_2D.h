//-*-C++-*-

/*
 *      Author: peter staar
 */


/*
#ifndef tight_binding_models_cluster_4A_2D_H_
#define tight_binding_models_cluster_4A_2D_H_

//#include "tight_binding_models.h"

//#include "No_symmetry.h"
//#include "D4.h"

template<>
class tight_binding_models<cluster_4A_2D> {

 public:

  // taken care off via parameters !
  static int               get_DCA_size();
  static int               get_LDA_size();

  static std::vector<int>& DCA_grid_size();
  static std::vector<int>& LDA_grid_size();

  // taken care of via model-type
  typedef Null_symmetry_2D  LDA_point_group;
  typedef D4                DCA_point_group;

  const static int DIMENSION = 2;
  static int       BANDS; //    = 1;

  const static double Fermi_energy = 0.;
  const static double density      = 1.;

  static double* get_r_DCA_basis();
  static double* get_k_DCA_basis();

  static double* get_r_LDA_basis();
  static double* get_k_LDA_basis();

  static std::complex<double> get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2);

  static void                 symmetrize_Hamiltonian(std::complex<double>* H_matrix);

  static std::vector<int>     get_interacting_bands();

  // taken care off via interaction-model-type
  static double               interaction_Hamiltonian(int nu_i, int nu_j, int delta_r=0);
  static double               alpha_field            (int nu_i, int nu_j, int delta_r=0);

  template<class parameters_type>
  static void   initialize_interaction_Hamiltonian(parameters_type& parameters);



private:
  static double U;
};

double tight_binding_models<cluster_4A_2D>::U = 0;
int    tight_binding_models<cluster_4A_2D>::BANDS = 1;

int tight_binding_models<cluster_4A_2D>::get_DCA_size()
{
  return DCA_grid_size()[0]*DCA_grid_size()[1];
}

int tight_binding_models<cluster_4A_2D>::get_LDA_size()
{
  return LDA_grid_size()[0]*LDA_grid_size()[1];
}

std::vector<int>& tight_binding_models<cluster_4A_2D>::DCA_grid_size()
{
  static std::vector<int> v(2);
  v[0] = 2;
  v[1] = 2;

  return v;
}


std::vector<int>& tight_binding_models<cluster_4A_2D>::LDA_grid_size()
{
  static std::vector<int> v(2);
  v[0] = 20;
  v[1] = 20;

  return v;
}


double* tight_binding_models<cluster_4A_2D>::get_r_DCA_basis()
{
  static double* r_DCA = new double[4];
  
  r_DCA[0] = 1.;  r_DCA[1] = 0.;
  r_DCA[2] = 0.;  r_DCA[3] = 1.;

  return r_DCA;
}


double* tight_binding_models<cluster_4A_2D>::get_k_DCA_basis()
{
   static double* k_DCA = new double[4];
  
  k_DCA[0] = 2*M_PI;  k_DCA[1] = 0.;
  k_DCA[2] = 0.;      k_DCA[3] = 2*M_PI;

  return k_DCA;
}


double* tight_binding_models<cluster_4A_2D>::get_r_LDA_basis()
{
  static double* r_LDA = new double[4];
  
  r_LDA[0] = 1.;  r_LDA[1] = 0.;
  r_LDA[2] = 0.;  r_LDA[3] = 1.;

  return r_LDA;
}


double* tight_binding_models<cluster_4A_2D>::get_k_LDA_basis()
{
  static double* k_LDA = new double[4];
  
  k_LDA[0] = 2.*M_PI;  k_LDA[1] = 0.;
  k_LDA[2] = 0.;       k_LDA[3] = 2.*M_PI;

  return k_LDA;
}


std::complex<double> tight_binding_models<cluster_4A_2D>::get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2)
{
  std::complex<double> H_LDA = 0.;

  if( (b1 == 0 && b2 == 0) && ((s1 == 0 && s2 == 0) || (s1 == 1 && s2 == 1)))
    H_LDA = -2.* t * ( cos(k[0]) + cos(k[1]) );

  return H_LDA;
}

void tight_binding_models<cluster_4A_2D>::symmetrize_Hamiltonian(std::complex<double>* H_matrix)
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

double tight_binding_models<cluster_4A_2D>::interaction_Hamiltonian(int nu_i, int nu_j, int delta_r)
{
  if( (nu_i==nu_j) && (delta_r==0) )
    return 0.;
  else
    return U;//U[delta_r];
}


template<class parameters_type>
void tight_binding_models<cluster_4A_2D>::initialize_interaction_Hamiltonian(parameters_type& parameters)
{
  U = parameters.get_U_hubbard();
}


std::vector<int> tight_binding_models<cluster_4A_2D>::get_interacting_bands()
{
  std::vector<int> v;
  v.push_back(0);
  return v;
}

#endif
*/
