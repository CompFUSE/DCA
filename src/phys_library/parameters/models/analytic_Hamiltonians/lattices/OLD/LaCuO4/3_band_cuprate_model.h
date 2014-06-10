//-*-C++-*-

/*
 *      Author: peter staar
 */

/*
#ifndef CUPRATE_3_BANDS_HAMILTONIANS
#define CUPRATE_3_BANDS_HAMILTONIANS

template<typename DCA_point_group_type>
class cuprate_model_3_bands
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
  static void initialize_H_interaction(function<double , domain >& H_interaction,
				       parameters_type&            parameters);

  template<class domain>
  static void initialize_H_symmetry(function<int , domain>& H_symmetry);

  template<class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters, 
						   std::vector<double> k, int b1, int s1, int b2, int s2);
};


template<typename DCA_point_group_type>
double* cuprate_model_3_bands<DCA_point_group_type>::initialize_r_DCA_basis()
{
  static double* r_DCA = new double[4];
  
  r_DCA[0] = 1.;  r_DCA[1] = 0.;
  r_DCA[2] = 0.;  r_DCA[3] = 1.;

  return r_DCA;
}

template<typename DCA_point_group_type>
double* cuprate_model_3_bands<DCA_point_group_type>::initialize_k_DCA_basis()
{
  static double* k_DCA = new double[4];
  
  k_DCA[0] = 2*M_PI;  k_DCA[1] = 0.;
  k_DCA[2] = 0.;      k_DCA[3] = 2*M_PI;

  return k_DCA;
}

template<typename DCA_point_group_type>
double* cuprate_model_3_bands<DCA_point_group_type>::initialize_r_LDA_basis()
{
  static double* r_LDA = new double[4];
  
  r_LDA[0] = 1.;  r_LDA[1] = 0.;
  r_LDA[2] = 0.;  r_LDA[3] = 1.;

  return r_LDA;
}

template<typename DCA_point_group_type>
double* cuprate_model_3_bands<DCA_point_group_type>::initialize_k_LDA_basis()
{
  static double* k_LDA = new double[4];
  
  k_LDA[0] = 2.*M_PI;  k_LDA[1] = 0.;
  k_LDA[2] = 0.;       k_LDA[3] = 2.*M_PI;

  return k_LDA;
}

template<typename DCA_point_group_type>
template<class domain, class parameters_type>
void cuprate_model_3_bands<DCA_point_group_type>::initialize_H_interaction(function<double , domain >& H_interaction,
									   parameters_type&            parameters)
{
  double U = parameters.get_U_hubbard();

  H_interaction(0,0,0,1,0) = U;
  H_interaction(0,1,0,0,0) = U;
}

template<typename DCA_point_group_type>
template<class domain>
void cuprate_model_3_bands<DCA_point_group_type>::initialize_H_symmetry(function<int , domain>& H_symmetries)
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
std::complex<double> cuprate_model_3_bands<DCA_point_group_type>::get_LDA_Hamiltonians(parameters_type& parameters,
										       std::vector<double> k, 
										       int b1, int s1, int b2, int s2)
{
  std::complex<double> i(0.,1.);
  std::complex<double> H_LDA = 0;

  if(s1==s2)
    {
      if(b1==0 && b2==0){
	H_LDA = parameters.get_dd()[0] 
	  + 2.*parameters.get_dd()[1]*(cos(k[0])+cos(k[1])) 
	  + 4.*parameters.get_dd()[2]*(cos(k[0])*cos(k[1]));
      }

      if (b1==1 && b2==1){
	H_LDA = (2.*parameters.get_pp()[1]*cos(k[0])
		 + 4.*parameters.get_pp()[2]*cos(k[0])*cos(k[1]));
      }

      if (b1==2 && b2==2){
	H_LDA = (2.*parameters.get_pp()[1]*cos(k[1])
		 + 4.*parameters.get_pd()[2]*cos(k[0])*cos(k[1]));
      }

      if ((b1==1 && b2==0) || (b1==0 && b2==1)){
	std::complex<double> result = (2.*parameters.get_pd()[0] + 4.*parameters.get_pd()[1]*cos(k[1])) * sin(k[0]/2.) * (-i) * exp(i*k[0]/2.);
	if(b1<b2)
	  H_LDA=result;
	else
	  H_LDA=conj(result);
      }

      if ((b1==2 && b2==0) || (b1==0 && b2==2)){
	std::complex<double> result = -(2.*parameters.get_pd()[0] + 4.*parameters.get_pd()[1]*cos(k[0])) * sin(k[1]/2.) * (-i) * exp(i*k[1]/2.);
	if(b1<b2)
	  H_LDA=result;
	else
	  H_LDA=conj(result);
      }

      if ((b1==1 && b2==2) || (b1==2 && b2==1)){
	std::complex<double> result = -( 4.*parameters.get_pp()[0]*sin(k[0]/2.)*sin(k[1]/2.))*exp(i*(k[1]-k[0])/2.);
	  
	if(b1<b2)
	  H_LDA=result;
	else
	  H_LDA=conj(result);
      }

    }
  
  return H_LDA;
}


#endif
*/
