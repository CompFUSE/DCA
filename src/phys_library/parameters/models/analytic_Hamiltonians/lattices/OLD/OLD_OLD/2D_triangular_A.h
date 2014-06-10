//-*-C++-*-

/*
 *      Author: peter staar
 */


#ifndef TRIANGULAR_A_2D_H_
#define TRIANGULAR_A_2D_H_



class triangular_A_2D
{
public:

  typedef Null_symmetry_2D LDA_point_group ;
  typedef /*Null_symmetry_2D*/ D6 DCA_point_group;

  const static int DIMENSION = 2;
  const static int BANDS     = 1;

  static double* initialize_r_DCA_basis();
  static double* initialize_k_DCA_basis();

  static double* initialize_r_LDA_basis();
  static double* initialize_k_LDA_basis();

  static std::complex<double> get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2);

  static void                 symmetrize_Hamiltonian(std::complex<double>* H_matrix);

  static std::vector<int>     initialize_interacting_bands();

  static std::vector<int>     get_interacting_bands();

};

double* triangular_A_2D::initialize_r_DCA_basis()
{
  static double* r_DCA = new double[4];
  
  r_DCA[0] = cos( M_PI/3.);  r_DCA[1] = sin( M_PI/3.);
  r_DCA[2] = cos(-M_PI/3.);  r_DCA[3] = sin(-M_PI/3.);

  return r_DCA;
}


double* triangular_A_2D::initialize_k_DCA_basis()
{
  static double* k_DCA = new double[4];
  
  k_DCA[0] = 2*M_PI;  k_DCA[1] = 2*M_PI/sqrt(3.);
  k_DCA[2] = 2*M_PI;  k_DCA[3] = -2*M_PI/sqrt(3.);

  return k_DCA;
}


double* triangular_A_2D::initialize_r_LDA_basis()
{
  static double* r_LDA = new double[4];
  
  r_LDA[0] = cos( M_PI/3.);  r_LDA[1] = sin( M_PI/3.);
  r_LDA[2] = cos(-M_PI/3.);  r_LDA[3] = sin(-M_PI/3.);

  return r_LDA;
}


double* triangular_A_2D::initialize_k_LDA_basis()
{
  static double* k_LDA = new double[4];
  
  k_LDA[0] = 2*M_PI;  k_LDA[1] = 2*M_PI/sqrt(3.);
  k_LDA[2] = 2*M_PI;  k_LDA[3] = -2*M_PI/sqrt(3.);

  return k_LDA;
}


std::complex<double> triangular_A_2D::get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2)
{
  std::complex<double> H_LDA = 0.;

  if( (b1 == 0 && b2 == 0) && ((s1 == 0 && s2 == 0) || (s1 == 1 && s2 == 1)))
    H_LDA = -2.* t * ( cos(k[0]) 
		       + cos(cos( M_PI/3.)*k[0] + sin( M_PI/3.)*k[1]) 
		       + cos(cos(-M_PI/3.)*k[0] + sin(-M_PI/3.)*k[1]));

  return H_LDA;
}

void triangular_A_2D::symmetrize_Hamiltonian(std::complex<double>* H_matrix)
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

std::vector<int> triangular_A_2D::get_interacting_bands()
{
  std::vector<int> v(1,0);
  return v;
}

#endif
