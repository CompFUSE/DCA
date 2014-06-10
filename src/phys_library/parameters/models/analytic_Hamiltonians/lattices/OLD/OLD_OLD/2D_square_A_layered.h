//-*-C++-*-

/*
 *      Author: peter staar
 */


#ifndef SQUARE_A_2D_LAYERED_H_
#define SQUARE_A_2D_LAYERED_H_



class square_A_2D_layered
{
public:

  typedef Null_symmetry_2D LDA_point_group ;
  typedef D4               DCA_point_group;

  const static int DIMENSION = 2;
  const static int BANDS     = 2;

  static double* initialize_r_DCA_basis();
  static double* initialize_k_DCA_basis();

  static double* initialize_r_LDA_basis();
  static double* initialize_k_LDA_basis();

  template<class domain, class parameters_type>
  static void initialize_H_interaction(function<double , domain >& H_interaction,
				       parameters_type&            parameters);

  static std::complex<double> get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2);

  static void                 symmetrize_Hamiltonian(std::complex<double>* H_matrix);

  static std::vector<int>     initialize_interacting_bands();

  static std::vector<int>     get_interacting_bands();
};

double* square_A_2D_layered::initialize_r_DCA_basis()
{
  static double* r_DCA = new double[4];
  
  r_DCA[0] = 1.;  r_DCA[1] = 0.;
  r_DCA[2] = 0.;  r_DCA[3] = 1.;

  return r_DCA;
}

double* square_A_2D_layered::initialize_k_DCA_basis()
{
   static double* k_DCA = new double[4];
  
  k_DCA[0] = 2*M_PI;  k_DCA[1] = 0.;
  k_DCA[2] = 0.;      k_DCA[3] = 2*M_PI;

  return k_DCA;
}

double* square_A_2D_layered::initialize_r_LDA_basis()
{
  static double* r_LDA = new double[4];
  
  r_LDA[0] = 1.;  r_LDA[1] = 0.;
  r_LDA[2] = 0.;  r_LDA[3] = 1.;

  return r_LDA;
}


double* square_A_2D_layered::initialize_k_LDA_basis()
{
  static double* k_LDA = new double[4];
  
  k_LDA[0] = 2.*M_PI;  k_LDA[1] = 0.;
  k_LDA[2] = 0.;       k_LDA[3] = 2.*M_PI;

  return k_LDA;
}

template<class domain, class parameters_type>
void square_A_2D_layered::initialize_H_interaction(function<double , domain >& H_interaction,
						   parameters_type&            parameters)
{
  double U = parameters.get_U_hubbard();
  double J = parameters.get_J_exchange();

  H_interaction(0,0,0) = 0;  H_interaction(0,1,0) = J;  H_interaction(0,2,0) = U;  H_interaction(0,3,0) = J;
  H_interaction(1,0,0) = J;  H_interaction(1,1,0) = 0;  H_interaction(1,2,0) = J;  H_interaction(1,3,0) = U;
  H_interaction(2,0,0) = U;  H_interaction(2,1,0) = J;  H_interaction(2,2,0) = 0;  H_interaction(2,3,0) = J;
  H_interaction(3,0,0) = J;  H_interaction(3,1,0) = U;  H_interaction(3,2,0) = J;  H_interaction(3,3,0) = 0;
}

std::complex<double> square_A_2D_layered::get_LDA_Hamiltonians(double t, std::vector<double> k, int b1, int s1, int b2, int s2)
{
  std::complex<double> H_LDA = 0.;

  if( (s1 == 0 && s2 == 0) || (s1 == 1 && s2 == 1) )
    {
      if( (b1 == b2) )
	H_LDA = -2.* t * ( cos(k[0]) + cos(k[1]) );
      else
	H_LDA = -0.1;
    }

  return H_LDA;
}

void square_A_2D_layered::symmetrize_Hamiltonian(std::complex<double>* H_matrix)
{
  //const static int matrix_size = (2*BANDS)*(2*BANDS);
  const static int matrix_dim  = (2*BANDS);

  // symmetrize the e_spin up and dn part
  std::complex<double> tmp_0_0 = (H_matrix[0+matrix_dim*0] + H_matrix[(2+0)+matrix_dim*(2+0)])/double(2); 
  std::complex<double> tmp_1_0 = (H_matrix[1+matrix_dim*0] + H_matrix[(2+1)+matrix_dim*(2+0)])/double(2); 
  std::complex<double> tmp_0_1 = (H_matrix[0+matrix_dim*1] + H_matrix[(2+0)+matrix_dim*(2+1)])/double(2); 
  std::complex<double> tmp_1_1 = (H_matrix[1+matrix_dim*1] + H_matrix[(2+1)+matrix_dim*(2+1)])/double(2); 

  // symmetry A <--> B
  H_matrix[0+matrix_dim*0] = (tmp_0_0 + tmp_1_1)/2.;
  H_matrix[1+matrix_dim*0] = (tmp_1_0 + tmp_0_1)/2.;
  H_matrix[0+matrix_dim*1] = (tmp_1_0 + tmp_0_1)/2.;
  H_matrix[1+matrix_dim*1] = (tmp_0_0 + tmp_1_1)/2.;

  H_matrix[(2+0)+matrix_dim*(2+0)] = (tmp_0_0 + tmp_1_1)/2.;
  H_matrix[(2+1)+matrix_dim*(2+0)] = (tmp_1_0 + tmp_0_1)/2.;
  H_matrix[(2+0)+matrix_dim*(2+1)] = (tmp_1_0 + tmp_0_1)/2.;
  H_matrix[(2+1)+matrix_dim*(2+1)] = (tmp_0_0 + tmp_1_1)/2.;
}

std::vector<int> square_A_2D_layered::get_interacting_bands()
{
  std::vector<int> v(0,0);
  v.push_back(0);
  v.push_back(1);
  return v;
}

#endif
