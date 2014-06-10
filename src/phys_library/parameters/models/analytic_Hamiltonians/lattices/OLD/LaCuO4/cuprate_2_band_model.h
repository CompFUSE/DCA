//-*-C++-*-

/*
 *      Author: peter staar
 */

#ifndef CUPRATE_2_BAND_MODEL_H
#define CUPRATE_2_BAND_MODEL_H

template<typename DCA_point_group_type>
class cuprate_2_band_model
{
public:

  //typedef Null_symmetry_2D     LDA_point_group;

  typedef no_symmetry<2>       LDA_point_group;
  typedef DCA_point_group_type DCA_point_group; 

  const static cluster_shape_type DCA_cluster_shape = BETT_CLUSTER;
  const static cluster_shape_type LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS     = 2;

  static double* initialize_r_DCA_basis();
  static double* initialize_k_DCA_basis();

  static double* initialize_r_LDA_basis();
  static double* initialize_k_LDA_basis();

  static std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > get_orbital_permutations();

  template<class domain, class parameters_type>
  static void initialize_H_interaction(function<double , domain >& H_interaction,
				       parameters_type&            parameters);

  template<class domain>
  static void initialize_H_symmetry(function<int , domain>& H_symmetry);

  template<class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters, 
						   std::vector<double> k, int b1, int s1, int b2, int s2);

private:

  static std::complex<double> beta(std::vector<double> k);
};


template<typename DCA_point_group_type>
double* cuprate_2_band_model<DCA_point_group_type>::initialize_r_DCA_basis()
{
  static double* r_DCA = new double[4];
  
  r_DCA[0] = 1.0;  r_DCA[1] = 0.0;
  r_DCA[2] = 0.0;  r_DCA[3] = 1.0;

  return r_DCA;
}

template<typename DCA_point_group_type>
double* cuprate_2_band_model<DCA_point_group_type>::initialize_k_DCA_basis()
{
  static double* k_DCA = new double[4];

  k_DCA[0] = 6.283184;  k_DCA[1] = 0.0;
  k_DCA[2] = 0.0;  k_DCA[3] = 6.283184;

  return k_DCA;
}

template<typename DCA_point_group_type>
double* cuprate_2_band_model<DCA_point_group_type>::initialize_r_LDA_basis()
{
  static double* r_LDA = new double[4];

  r_LDA[0] = 0.5;  r_LDA[1] = 0.0;
  r_LDA[2] = 0.0;  r_LDA[3] = 0.5;

  return r_LDA;
}

template<typename DCA_point_group_type>
double* cuprate_2_band_model<DCA_point_group_type>::initialize_k_LDA_basis()
{
  static double* k_LDA = new double[4];
  
  k_LDA[0] = 4*M_PI;  k_LDA[1] = 0.0;
  k_LDA[2] = 0.0;     k_LDA[3] = 4*M_PI;

  return k_LDA;
}

template<typename DCA_point_group_type>
template<class domain, class parameters_type>
void cuprate_2_band_model<DCA_point_group_type>::initialize_H_interaction(function<double , domain >& H_interaction,
						      	  parameters_type&            parameters)
{
  std::vector<std::vector<double> >& U_ij = parameters.get_U_ij();

  for(size_t i=0; i<U_ij.size(); i++)
    H_interaction(U_ij[i][0], U_ij[i][1], U_ij[i][2], U_ij[i][3], U_ij[i][4]) = U_ij[i][5];
}

template<typename DCA_point_group_type>
std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > cuprate_2_band_model<DCA_point_group_type>::get_orbital_permutations()
{
  std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > permutations(0);
  return permutations;
}

template<typename DCA_point_group_type>
template<class domain>
void cuprate_2_band_model<DCA_point_group_type>::initialize_H_symmetry(function<int , domain>& H_symmetries)
{
  H_symmetries = -1.;

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
std::complex<double> cuprate_2_band_model<DCA_point_group_type>::get_LDA_Hamiltonians(parameters_type& parameters,
										      std::vector<double> k, 
										      int b1, int s1, int b2, int s2)
{
  std::vector<std::vector<double> >& t_ij = parameters.get_t_ij();

//   double t=0;
//   for(size_t i=0; i<t_ij.size(); i++)
//     if(t_ij[i][0]==b1 && t_ij[i][1]==b2 && t_ij[i][2]==0)
//       t = t_ij[i][3];

  static double E_d, E_p, t_pp, t_pd;

  for(size_t i=0; i<t_ij.size(); i++)
    {
      if(t_ij[i][0]==0 && t_ij[i][1]==0 && t_ij[i][2]==0)
	E_d = t_ij[i][3];
      
      if(t_ij[i][0]==1 && t_ij[i][1]==1 && t_ij[i][2]==0)
	E_p = t_ij[i][3];
      
      if(((t_ij[i][0]==1 && t_ij[i][1]==0) || (t_ij[i][0]==0 && t_ij[i][1]==1)) && t_ij[i][2]==0)
	t_pd = t_ij[i][3];
      
      if(((t_ij[i][0]==1 && t_ij[i][1]==0) || (t_ij[i][0]==0 && t_ij[i][1]==1)) && t_ij[i][2]==1)
	t_pp = t_ij[i][3];

//       cout << E_d << "\t" << E_p << "\t" << t_pd << "\t" << t_pp << "\n"; 
    }

  std::complex<double> H_LDA = 0;

  if(s1==s2)
  {
    if(b1==0 && b2==0)
      H_LDA = E_d;

    if((b1==0 && b2==1) || (b1==1 && b2==0))
      H_LDA = 2.*t_pd*sqrt(square(sin(k[0]/2.)) + square(sin(k[1]/2.)));

    if(b1==1 && b2==1)
      H_LDA = E_p - 8.*t_pp*square(beta(k)*sin(k[0]/2.)*sin(k[1]/2.));

//     cout << H_LDA << endl;
 }

  return H_LDA;
}

template<typename DCA_point_group_type>
std::complex<double> cuprate_2_band_model<DCA_point_group_type>::beta(std::vector<double> k)
{
  if(fabs(k[0])<1.e-3 && fabs(k[1])<1.e-3)
    return 0.;
  else
    return 1./sqrt(square(sin(k[0]/2.)) + square(sin(k[1]/2.)));
}

#endif
