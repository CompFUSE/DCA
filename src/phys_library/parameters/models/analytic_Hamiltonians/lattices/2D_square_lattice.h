//-*-C++-*-

#ifndef SQUARE_LATTICE_2D_H
#define SQUARE_LATTICE_2D_H

/*!
 *  \author peter staar
 */
template<typename point_group_type>
class square_lattice
{
public:

  typedef no_symmetry<2>   LDA_point_group;
  typedef point_group_type DCA_point_group;

//   const static cluster_shape_type DCA_cluster_shape = BETT_CLUSTER;
//   const static cluster_shape_type LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS     = 1;

  static double* initialize_r_DCA_basis();
//   static double* initialize_k_DCA_basis();

  static double* initialize_r_LDA_basis();
//   static double* initialize_k_LDA_basis();

  static std::vector<int>                  get_flavors();
  static std::vector<std::vector<double> > get_a_vectors();

  static std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > get_orbital_permutations();

  template<class domain, class parameters_type>
  static void initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
                                       parameters_type&            parameters);

  template<class domain>
  static void initialize_H_symmetry(FUNC_LIB::function<int , domain>& H_symmetry);

  template<class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters,
                                                   std::vector<double> k, int b1, int s1, int b2, int s2);
};

template<typename point_group_type>
double* square_lattice<point_group_type>::initialize_r_DCA_basis()
{
  static double* r_DCA = new double[4];

  r_DCA[0] = 1.;  r_DCA[1] = 0.;
  r_DCA[2] = 0.;  r_DCA[3] = 1.;

  return r_DCA;
}

// template<typename point_group_type>
// double* square_lattice<point_group_type>::initialize_k_DCA_basis()
// {
//   static double* k_DCA = new double[4];

//   k_DCA[0] = 2*M_PI;  k_DCA[1] = 0.;
//   k_DCA[2] = 0.;      k_DCA[3] = 2*M_PI;

//   return k_DCA;
// }

template<typename point_group_type>
double* square_lattice<point_group_type>::initialize_r_LDA_basis()
{
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;  r_LDA[1] = 0.;
  r_LDA[2] = 0.;  r_LDA[3] = 1.;

  return r_LDA;
}

// template<typename point_group_type>
// double* square_lattice<point_group_type>::initialize_k_LDA_basis()
// {
//   static double* k_LDA = new double[4];

//   k_LDA[0] = 2.*M_PI;  k_LDA[1] = 0.;
//   k_LDA[2] = 0.;       k_LDA[3] = 2.*M_PI;

//   return k_LDA;
// }

template<typename point_group_type>
std::vector<int> square_lattice<point_group_type>::get_flavors()
{
  static std::vector<int> flavors(BANDS);

  for(int i=0; i<BANDS; i++)
    flavors[i]=i;

  return flavors;
}

template<typename point_group_type>
std::vector<std::vector<double> > square_lattice<point_group_type>::get_a_vectors()
{
  static std::vector<std::vector<double> > a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));
  return a_vecs;
}

template<typename point_group_type>
std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > square_lattice<point_group_type>::get_orbital_permutations()
{
  static std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > permutations(0);
  return permutations;
}

// TODO: Add non-local interaction of same spins.
//       Use V instead of U_prime?
template<typename point_group_type>
template<class domain, class parameters_type>
void square_lattice<point_group_type>::initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
                                                                parameters_type&            parameters)
{
  H_interaction = 0.;

  double U = parameters.get_U();
  double U_prime = parameters.get_U_prime();

  // actually the same as DCA_r_cluster_type (see typedifinitions.h).
  typedef typename TypeAt<typename domain::domain_typelist_2, 0>::Result DCA_r_cluster_t;

  int DIMENSION = DCA_r_cluster_t::DIMENSION;
  assert(DIMENSION == 2);
  
  int origin = DCA_r_cluster_t::origin_index();

  std::vector<typename DCA_r_cluster_t::element_type>& basis = DCA_r_cluster_t::get_basis_vectors();
  std::vector<typename DCA_r_cluster_t::element_type>& super_basis = DCA_r_cluster_t::get_super_basis_vectors();
  std::vector<typename DCA_r_cluster_t::element_type>& elements = DCA_r_cluster_t::get_elements();

  std::vector<int> nn_index(DIMENSION);  // Indices of nearest neighbours w.r.t. origin.
  for(int d = 0; d < DIMENSION; ++d) {
    std::vector<double> basis_vec = cluster_operations::translate_inside_cluster(basis[d], super_basis);
    nn_index[d] = cluster_operations::index(basis_vec, elements, BRILLOUIN_ZONE);
  }

  // non-local opposite-spin interaction
  H_interaction(0, 1, nn_index[0]) = U_prime;
  H_interaction(1, 0, nn_index[0]) = U_prime;

  H_interaction(0, 1, nn_index[1]) = U_prime;
  H_interaction(1, 0, nn_index[1]) = U_prime;

  // non-local same-spin interaction
  H_interaction(0, 0, nn_index[0]) = U_prime;
  H_interaction(1, 1, nn_index[0]) = U_prime;

  H_interaction(0, 0, nn_index[1]) = U_prime;
  H_interaction(1, 1, nn_index[1]) = U_prime;

  // local interaction
  H_interaction(0, 1, origin) = U;
  H_interaction(1, 0, origin) = U;
}

template<typename point_group_type>
template<class domain>
void square_lattice<point_group_type>::initialize_H_symmetry(FUNC_LIB::function<int , domain>& H_symmetries)
{
  H_symmetries(0,0)= 0; H_symmetries(0,1)=-1;
  H_symmetries(1,0)=-1; H_symmetries(1,1)= 0;
}

template<typename point_group_type>
template<class parameters_type>
std::complex<double> square_lattice<point_group_type>::get_LDA_Hamiltonians(parameters_type& parameters,
                                                                            std::vector<double> k,
                                                                            int b1, int s1, int b2, int s2)
{
  std::complex<double> H_LDA = 0.;

  double t       = parameters.get_t();
  double t_prime = parameters.get_t_prime();

  if( (b1 == 0 && b2 == 0) && ((s1 == 0 && s2 == 0) || (s1 == 1 && s2 == 1)))
    H_LDA = -2.* t * ( cos(k[0]) + cos(k[1]) ) - 4.*t_prime*cos(k[0])*cos(k[1]);

  return H_LDA;
}

#endif
