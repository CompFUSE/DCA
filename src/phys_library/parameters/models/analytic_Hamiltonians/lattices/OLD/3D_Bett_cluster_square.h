//-*-C++-*-

/*
 *      Author: peter staar
 */


#ifndef BETT_CLUSTER_SQUARE_3D
#define BETT_CLUSTER_SQUARE_3D

template<typename DCA_point_group_type>
class Bett_cluster_square_3D
{
public:

  //typedef Null_symmetry_2D     LDA_point_group;

  typedef no_symmetry<2>       LDA_point_group;
  typedef DCA_point_group_type DCA_point_group; 

  const static cluster_shape_type DCA_cluster_shape = BETT_CLUSTER;
  const static cluster_shape_type LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 3;
  const static int BANDS     = 1;

  static double* initialize_r_DCA_basis();
  static double* initialize_k_DCA_basis();

  static double* initialize_r_LDA_basis();
  static double* initialize_k_LDA_basis();

  static std::vector<int>                  get_flavors();
  static std::vector<std::vector<double> > get_a_vectors();

  template<class domain, class parameters_type>
  static void initialize_H_interaction(function<double , domain >& H_interaction,
				       parameters_type&            parameters);

  template<class domain>
  static void initialize_H_symmetry(function<int , domain>& H_symmetry);

  template<typename parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters, std::vector<double> k, int b1, int s1, int b2, int s2);

//   static void                 symmetrize_Hamiltonian(std::complex<double>* H_matrix);

//   static std::vector<int>     initialize_interacting_bands();

//   static std::vector<int>     get_interacting_bands();
};

template<typename DCA_point_group_type>
double* Bett_cluster_square_3D<DCA_point_group_type>::initialize_r_DCA_basis()
{
  static double* r_DCA = new double[9];
  
  r_DCA[0] = 1.;  r_DCA[1] = 0.;  r_DCA[2] = 0.;
  r_DCA[3] = 0.;  r_DCA[4] = 1.;  r_DCA[5] = 0.;
  r_DCA[6] = 0.;  r_DCA[7] = 0.;  r_DCA[8] = 1.;

  return r_DCA;
}

template<typename DCA_point_group_type>
double* Bett_cluster_square_3D<DCA_point_group_type>::initialize_k_DCA_basis()
{
  static double* k_DCA = new double[9];

  k_DCA[0] = 2*M_PI;  k_DCA[1] = 0.;      k_DCA[2] = 0.;
  k_DCA[3] = 0.;      k_DCA[4] = 2*M_PI;  k_DCA[5] = 0.;
  k_DCA[6] = 0.;      k_DCA[7] = 0.;      k_DCA[8] = 2*M_PI;

  return k_DCA;
}

template<typename DCA_point_group_type>
double* Bett_cluster_square_3D<DCA_point_group_type>::initialize_r_LDA_basis()
{
  static double* r_LDA = new double[9];
  
  r_LDA[0] = 1.;  r_LDA[1] = 0.;  r_LDA[2] = 0.;
  r_LDA[3] = 0.;  r_LDA[4] = 1.;  r_LDA[5] = 0.;
  r_LDA[6] = 0.;  r_LDA[7] = 0.;  r_LDA[8] = 1.;

  return r_LDA;
}

template<typename DCA_point_group_type>
double* Bett_cluster_square_3D<DCA_point_group_type>::initialize_k_LDA_basis()
{
  static double* k_LDA = new double[9];

  k_LDA[0] = 2*M_PI;  k_LDA[1] = 0.;      k_LDA[2] = 0.;
  k_LDA[3] = 0.;      k_LDA[4] = 2*M_PI;  k_LDA[5] = 0.;
  k_LDA[6] = 0.;      k_LDA[7] = 0.;      k_LDA[8] = 2*M_PI;

  return k_LDA;
}

template<typename DCA_point_group_type>
std::vector<int> Bett_cluster_square_3D<DCA_point_group_type>::get_flavors()
{
  static std::vector<int> flavors(BANDS);

  for(int i=0; i<BANDS; i++)
    flavors[i]=i;

  return flavors;
}

template<typename DCA_point_group_type>
std::vector<std::vector<double> > Bett_cluster_square_3D<DCA_point_group_type>::get_a_vectors()
{
  static std::vector<std::vector<double> > a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));
  return a_vecs;
}

template<typename DCA_point_group_type>
template<class domain, class parameters_type>
void Bett_cluster_square_3D<DCA_point_group_type>::initialize_H_interaction(function<double , domain >& H_interaction,
					   parameters_type&            parameters)
{
  std::vector<std::vector<double> >& U_ij = parameters.get_U_ij();

  for(size_t i=0; i<U_ij.size(); i++)
    H_interaction(U_ij[i][0], U_ij[i][1], U_ij[i][2], U_ij[i][3], U_ij[i][4]) = U_ij[i][5];

//   double U = parameters.get_U_hubbard();

//   H_interaction(0,0,0) = 0;  H_interaction(0,1,0) = U;
//   H_interaction(1,0,0) = U;  H_interaction(1,1,0) = 0;
}

template<typename DCA_point_group_type>
template<class domain>
void Bett_cluster_square_3D<DCA_point_group_type>::initialize_H_symmetry(function<int , domain>& H_symmetries)
{
  H_symmetries(0,0)= 0; H_symmetries(0,1)=-1;
  H_symmetries(1,0)=-1; H_symmetries(1,1)= 0;
}

template<typename DCA_point_group_type>
template<typename parameters_type>
std::complex<double> Bett_cluster_square_3D<DCA_point_group_type>::get_LDA_Hamiltonians(parameters_type& parameters, std::vector<double> k, int b1, int s1, int b2, int s2)
{
  std::vector<std::vector<double> >& t_ij = parameters.get_t_ij();

  double t=0;
  for(size_t i=0; i<t_ij.size(); i++)
    if(t_ij[i][0]==b1 && t_ij[i][1]==b2 && t_ij[i][2]==0)
      t = t_ij[i][3];

  std::complex<double> H_LDA = 0.;

  if( (b1 == 0 && b2 == 0) && ((s1 == 0 && s2 == 0) || (s1 == 1 && s2 == 1)))
    H_LDA = -2.* t * ( cos(k[0]) + cos(k[1]) + cos(k[2]) );

  return H_LDA;
}

#endif
