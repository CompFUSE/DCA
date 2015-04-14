//-*-C++-*-

/*
 *      Author: peter staar
 */

#ifndef ANDERSEN_HAMILTONIANS
#define ANDERSEN_HAMILTONIANS

template<typename DCA_point_group_type>
class Andersen_Hamiltonians
{
public:

  typedef Null_symmetry_2D     LDA_point_group;
  typedef DCA_point_group_type DCA_point_group; 

  const static cluster_shape_type DCA_cluster_shape = BETT_CLUSTER;
  const static cluster_shape_type LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS     = 3;

  static double* initialize_r_DCA_basis();
  static double* initialize_k_DCA_basis();

  static double* initialize_r_LDA_basis();
  static double* initialize_k_LDA_basis();

  template<class domain, class parameters_type>
  static void initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
				       parameters_type&            parameters);

  template<class domain>
  static void initialize_H_symmetry(FUNC_LIB::function<int , domain>& H_symmetry);

  template<class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters, 
						   std::vector<double> k, int b1, int s1, int b2, int s2);

private:

  static std::complex<double> Edk(std::vector<double> k);
  static std::complex<double> E1k(std::vector<double> k);
  static std::complex<double> E2k(std::vector<double> k);

  static std::complex<double> V1dk(std::vector<double> k);
  static std::complex<double> V2dk(std::vector<double> k);
  static std::complex<double> V12k(std::vector<double> k);

private:

  static const double twor;
  static const double fourr;
  static const double halfr;

  static const double h_tdd1;
  static const double h_tdd2;

  static const double h_txx1;
  static const double h_txx2;
  static const double h_txx3;

  static const double h_txd1;
  static const double h_txd2;
  
  static const double h_txy1;

  // Hopping parameters for H(k)

  // Taken from "3-orbital LDA Hamiltonian" PDF by O. K. Anderson,
  // T. Saha-Dasgupta. O. Jepsen revision 20th February 2004

  // We follow the sign convention and paramter labels in the above
  // paper. Energy units are eV.

  // Example:

  // 8<-----La2CuO4 parameters-------
  // 0.431   ! h_edmep e_d-e_p                                        
  // -0.102  ! h_tdd1  d(000)->d(100)                                  
  // 0.034   ! h_tdd2  d(000)->d(110)                                  
  // 0.961   ! h_txd1  d(000)->x(1/2 0 0)                              
  // -0.103  ! h_txd2  d(000)->x(1/2 1 0)                              
  // -0.001  ! h_txd3  d(000)->x(3/2 0 0)                              
  // 0.022   ! h_txd4  d(000)->x(3/2 1 0)                              
  // 0.001   ! h_txd5  d(000)->x(1/2 2 0)                              
  // -0.005  ! h_txd6  d(000)->x(3/2 2 0)                              
  // 0.001   ! h_txd7  d(000)->x(5/2 0 0)                              
  // 0.154   ! h_txy1  x(-1/2 0 0)->y(0 1/2 0)                         
  // -0.026  ! h_txy2  x(-1/2 0 0)->y(1 1/2 0)                         
  // 0.018   ! h_txy3  x(-1/2 0 0)->y(1 3/2 0)                         
  // -0.004  ! h_txy4  x(-1/2 0 0)->y(1 5/2 0)                         
  // 0.034   ! h_txy5  x(-1/2 0 0)->y(1/2 1 1) ONLY in bct compounds   
  // 0.0     ! h_txy6  x(-1/2 0 0)->y(0 1/2 1) ONLY in sc compounds    
  // -0.241  ! h_txx1  x(-1/2 0 0)->x(1/2 0 0)                         
  // 0.019   ! h_txx2  x( 1/2 0 0)->x(1/2 1 0)                         
  // 0.107   ! h_txx3  x(-1/2 0 0)->x(1/2 1 0)                         
  // -0.019  ! h_txx4  x(-1/2 0 0)->x(3/2 0 0)                         
  // -0.020  ! h_txx5  x(-1/2 0 0)->x(1/2 2 0)                         
  // -0.010  ! h_txx6  x(-1/2 0 0)->x(3/2 1 0)                         
  // 0.003   ! h_txx7  x(-1/2 0 0)->x(3/2 2 0)                         
  // -0.019  ! h_txx8  x( 1/2 0 0)->x(1/2 2 0)                         
  // 0.007   ! h_txx9  x( 1/2 0 0)->x(1/2 3 0)                         
  // -0.001  ! h_txx10 x(-1/2 0 0)->x(5/2 0 0)                        
  // 8<-------------------------------
};

template<typename DCA_point_group_type>
const double Andersen_Hamiltonians<DCA_point_group_type>::twor  = 2.;

template<typename DCA_point_group_type>
const double Andersen_Hamiltonians<DCA_point_group_type>::fourr = 4.;

template<typename DCA_point_group_type>
const double Andersen_Hamiltonians<DCA_point_group_type>::halfr = 1./2.;

template<typename DCA_point_group_type>
const double Andersen_Hamiltonians<DCA_point_group_type>::h_tdd1 = -0.102;
template<typename DCA_point_group_type>
const double Andersen_Hamiltonians<DCA_point_group_type>::h_tdd2 =  0.034;

template<typename DCA_point_group_type>
const double Andersen_Hamiltonians<DCA_point_group_type>::h_txx1 = -0.241;
template<typename DCA_point_group_type>
const double Andersen_Hamiltonians<DCA_point_group_type>::h_txx2 =  0.019;
template<typename DCA_point_group_type>
const double Andersen_Hamiltonians<DCA_point_group_type>::h_txx3 =  0.107;

template<typename DCA_point_group_type>
const double Andersen_Hamiltonians<DCA_point_group_type>::h_txd1 =  0.961;
template<typename DCA_point_group_type>
const double Andersen_Hamiltonians<DCA_point_group_type>::h_txd2 = -0.103;
template<typename DCA_point_group_type>  
const double Andersen_Hamiltonians<DCA_point_group_type>::h_txy1 = 0.154;


template<typename DCA_point_group_type>
double* Andersen_Hamiltonians<DCA_point_group_type>::initialize_r_DCA_basis()
{
  static double* r_DCA = new double[4];
  
  r_DCA[0] = 1.;  r_DCA[1] = 0.;
  r_DCA[2] = 0.;  r_DCA[3] = 1.;

  return r_DCA;
}

template<typename DCA_point_group_type>
double* Andersen_Hamiltonians<DCA_point_group_type>::initialize_k_DCA_basis()
{
  static double* k_DCA = new double[4];
  
  k_DCA[0] = 2*M_PI;  k_DCA[1] = 0.;
  k_DCA[2] = 0.;      k_DCA[3] = 2*M_PI;

  return k_DCA;
}

template<typename DCA_point_group_type>
double* Andersen_Hamiltonians<DCA_point_group_type>::initialize_r_LDA_basis()
{
  static double* r_LDA = new double[4];
  
  r_LDA[0] = 1.;  r_LDA[1] = 0.;
  r_LDA[2] = 0.;  r_LDA[3] = 1.;

  return r_LDA;
}

template<typename DCA_point_group_type>
double* Andersen_Hamiltonians<DCA_point_group_type>::initialize_k_LDA_basis()
{
  static double* k_LDA = new double[4];
  
  k_LDA[0] = 2.*M_PI;  k_LDA[1] = 0.;
  k_LDA[2] = 0.;       k_LDA[3] = 2.*M_PI;

  return k_LDA;
}

template<typename DCA_point_group_type>
template<class domain, class parameters_type>
void Andersen_Hamiltonians<DCA_point_group_type>::initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
										    parameters_type&            parameters)
{
  double U = parameters.get_U_hubbard();

  H_interaction(0,0,0,1,0) = U;
  H_interaction(0,1,0,0,0) = U;
}

template<typename DCA_point_group_type>
template<class domain>
void Andersen_Hamiltonians<DCA_point_group_type>::initialize_H_symmetry(FUNC_LIB::function<int , domain>& H_symmetries)
{
  // e_up <==> e_dn symmetry
  for(int i=0; i<BANDS; i++){
    for(int j=0; j<BANDS; j++){
      int l = j+i*BANDS;
      H_symmetries(i,0,j,0) = l;
      H_symmetries(i,1,j,1) = l;
    }
  }
}

template<typename DCA_point_group_type>
template<class parameters_type>
std::complex<double> Andersen_Hamiltonians<DCA_point_group_type>::get_LDA_Hamiltonians(parameters_type& parameters,
												std::vector<double> k, 
												int b1, int s1, int b2, int s2)
{
  // Evaluate components of H(k) where 

  //       | Edk     V1dk   V2dk |
  //  H(k)=| V1dk*   E1k    V12k |
  //       | V2dk*   V12k*  E2k  |

  //  And H(k) is hermitian
  //  Hole representation
  //  Convention for hoppins determination |k>=sum_i e^(-ii k r_i)|i>
  //  Minus signs DO NOT follow O. K. Anderson

  std::complex<double> H_LDA = 0.;
  //double t       = parameters.get_t_hopping();
  //double t_prime = parameters.get_t_prime_hopping();

  if( ((s1 == 0 && s2 == 0) || (s1 == 1 && s2 == 1)))
    {
      if(b1==0 && b2==0)
	H_LDA = Edk(k);

      if(b1==1 && b2==1)
	H_LDA = E1k(k);

      if(b1==2 && b2==2)
	H_LDA = E2k(k);

      if(b1==1 && b2==0)
	H_LDA = conj(V1dk(k));

      if(b1==2 && b2==0)
	H_LDA = conj(V2dk(k));

      if(b1==2 && b2==1)
	H_LDA = conj(V12k(k));

      if(b1==0 && b2==1)
	H_LDA = V1dk(k);

      if(b1==0 && b2==2)
	H_LDA = V2dk(k);

      if(b1==1 && b2==2)
	H_LDA = V12k(k);
    }
  else
    H_LDA = 0.;

  return H_LDA;
}

template<typename DCA_point_group_type>
std::complex<double> Andersen_Hamiltonians<DCA_point_group_type>::Edk(std::vector<double> k)
{
  return twor*h_tdd1*(cos(k[0])+cos(k[1])) + fourr*h_tdd2*cos(k[0])*cos(k[1]);
}

template<typename DCA_point_group_type>
std::complex<double> Andersen_Hamiltonians<DCA_point_group_type>::E1k(std::vector<double> k)
{
  return twor*h_txx1*cos(k[0])+twor*h_txx2*cos(k[1]) + fourr*h_txx3*cos(k[0])*cos(k[1]);
}

template<typename DCA_point_group_type>
std::complex<double> Andersen_Hamiltonians<DCA_point_group_type>::E2k(std::vector<double> k)
{
  // Same dispersion as E1k, except kx<->ky
  return twor*h_txx1*cos(k[1])+twor*h_txx2*cos(k[0]) + fourr*h_txx3*cos(k[1])*cos(k[0]);
}

template<typename DCA_point_group_type>
std::complex<double> Andersen_Hamiltonians<DCA_point_group_type>::V1dk(std::vector<double> k)
{
  const static std::complex<double> ii(0.,1.);
  return   -twor*ii*sin(halfr*k[0])*h_txd1 - fourr*ii*sin(halfr*k[0])*cos(k[1])*h_txd2;
}

template<typename DCA_point_group_type>
std::complex<double> Andersen_Hamiltonians<DCA_point_group_type>::V2dk(std::vector<double> k)
{
  // Same as V1dk except kx<->ky and overall -ve
  const static std::complex<double> ii(0.,1.);
  return -(-twor*ii*sin(halfr*k[1])*h_txd1 - fourr*ii*sin(halfr*k[1])*cos(k[0])*h_txd2);
}

template<typename DCA_point_group_type>
std::complex<double> Andersen_Hamiltonians<DCA_point_group_type>::V12k(std::vector<double> k)
{
  return fourr*h_txy1*sin(halfr*k[0])*sin(halfr*k[1]);
}

#endif
