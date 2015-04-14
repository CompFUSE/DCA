//-*-C++-*-

/*
 *      Author: peter staar
 */


#ifndef SQUARE_A_3D_H_
#define SQUARE_A_3D_H_



class square_A_3D
{
public:

  typedef Null_symmetry_3D LDA_point_group ;
  typedef Oh               DCA_point_group;

  const static int DIMENSION = 3;
  const static int BANDS     = 1;

  static double* initialize_r_DCA_basis();
  static double* initialize_k_DCA_basis();

  static double* initialize_r_LDA_basis();
  static double* initialize_k_LDA_basis();

  template<class domain, class parameters_type>
  static void initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
				       parameters_type&            parameters);

  static std::complex<double> get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2);

  static void                 symmetrize_Hamiltonian(std::complex<double>* H_matrix);

  static std::vector<int>     initialize_interacting_bands();

  static std::vector<int>     get_interacting_bands();
};

double* square_A_3D::initialize_r_DCA_basis()
{
  static double* r_DCA = new double[9];
  
  r_DCA[0] = 1.;  r_DCA[1] = 0.;  r_DCA[2] = 0.;
  r_DCA[3] = 0.;  r_DCA[4] = 1.;  r_DCA[5] = 0.;
  r_DCA[6] = 0.;  r_DCA[7] = 0.;  r_DCA[8] = 1.;

  return r_DCA;
}

double* square_A_3D::initialize_k_DCA_basis()
{
  static double* k_DCA = new double[9];

  k_DCA[0] = 2*M_PI;  k_DCA[1] = 0.;      k_DCA[2] = 0.;
  k_DCA[3] = 0.;      k_DCA[4] = 2*M_PI;  k_DCA[5] = 0.;
  k_DCA[6] = 0.;      k_DCA[7] = 0.;      k_DCA[8] = 2*M_PI;

  return k_DCA;
}

double* square_A_3D::initialize_r_LDA_basis()
{
  static double* r_LDA = new double[9];
  
  r_LDA[0] = 1.;  r_LDA[1] = 0.;  r_LDA[2] = 0.;
  r_LDA[3] = 0.;  r_LDA[4] = 1.;  r_LDA[5] = 0.;
  r_LDA[6] = 0.;  r_LDA[7] = 0.;  r_LDA[8] = 1.;

  return r_LDA;
}

double* square_A_3D::initialize_k_LDA_basis()
{
  static double* k_LDA = new double[9];

  k_LDA[0] = 2*M_PI;  k_LDA[1] = 0.;      k_LDA[2] = 0.;
  k_LDA[3] = 0.;      k_LDA[4] = 2*M_PI;  k_LDA[5] = 0.;
  k_LDA[6] = 0.;      k_LDA[7] = 0.;      k_LDA[8] = 2*M_PI;

  return k_LDA;
}


template<class domain, class parameters_type>
void square_A_3D::initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
					   parameters_type&            parameters)
{
  double U = parameters.get_U_hubbard();

  H_interaction(0,0,0) = 0;  H_interaction(0,1,0) = U;
  H_interaction(1,0,0) = U;  H_interaction(1,1,0) = 0;
}


std::complex<double> square_A_3D::get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2)
{
  std::complex<double> H_LDA = 0.;

  if( (b1 == 0 && b2 == 0) && ((s1 == 0 && s2 == 0) || (s1 == 1 && s2 == 1)))
    H_LDA = -2.* t * ( cos(k[0]) + cos(k[1]) + cos(k[2]) );

  return H_LDA;
}

void square_A_3D::symmetrize_Hamiltonian(std::complex<double>* H_matrix)
{
  const static int matrix_dim  = (2*BANDS);

  // symmetrize the e_spin up and dn part
  std::complex<double> tmp = (H_matrix[0+matrix_dim*0] + H_matrix[1+matrix_dim*1])/double(2); 

  H_matrix[0+matrix_dim*0] = tmp;
  H_matrix[1+matrix_dim*0] = 0;
  H_matrix[0+matrix_dim*1] = 0;
  H_matrix[1+matrix_dim*1] = tmp;
}

std::vector<int> square_A_3D::get_interacting_bands()
{
  std::vector<int> v(1,0);
  return v;
}


#endif
