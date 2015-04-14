//-*-C++-*-

#ifndef BILAYER_MODEL_H
#define BILAYER_MODEL_H

/*!
 *  \author peter staar
 */
template<typename DCA_point_group_type>
class bilayer_model
{
  //#include type_definitions.h"
public:

  typedef no_symmetry<2>       LDA_point_group;
  typedef DCA_point_group_type DCA_point_group; 

  const static cluster_shape_type DCA_cluster_shape = BETT_CLUSTER;
  const static cluster_shape_type LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS     = 2;

public:

  template<class parameters_type>
  static void initialize(parameters_type& parameters);

  static std::vector<int>& LDA_grid_size();

  static double* get_r_DCA_basis();
  static double* get_k_DCA_basis();

  static double* get_r_LDA_basis();
  static double* get_k_LDA_basis();

  static double* initialize_r_DCA_basis();
  static double* initialize_k_DCA_basis();

  static double* initialize_r_LDA_basis();
  static double* initialize_k_LDA_basis();

  static std::vector<int>                  get_flavors();
  static std::vector<std::vector<double> > get_a_vectors();

  static std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > get_orbital_permutations();

  template<class domain, class parameters_type>
  static void initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
				       parameters_type&            parameters);

  template<class domain>
  static void initialize_H_symmetries(FUNC_LIB::function<int, domain>& H_symmetry);

  template<class domain, class parameters_type>
  static void initialize_H_LDA(FUNC_LIB::function<std::complex<double> , domain >& H_LDA,
			       parameters_type&                          parameters);

private:

  template<class parameters_type>
  static std::complex<double> get_LDA_Hamiltonians(parameters_type& parameters, 
						   std::vector<double> k, int b1, int s1, int b2, int s2);

  static double t_x;
  static double t_y;
  static double t_z;

  static double U;
  static double V;
};

template<typename DCA_point_group_type>
double bilayer_model<DCA_point_group_type>::t_x = 0;

template<typename DCA_point_group_type>
double bilayer_model<DCA_point_group_type>::t_y = 0;

template<typename DCA_point_group_type>
double bilayer_model<DCA_point_group_type>::t_z = 0;

template<typename DCA_point_group_type>
double bilayer_model<DCA_point_group_type>::U = 0;

template<typename DCA_point_group_type>
double bilayer_model<DCA_point_group_type>::V = 0;

template<typename DCA_point_group_type>
template<class parameters_type>
void bilayer_model<DCA_point_group_type>::initialize(parameters_type& parameters)
{
  t_x = parameters.get_tx();
  t_y = parameters.get_ty();
  t_z = parameters.get_tz();

  U = parameters.get_U();
  V = parameters.get_V();

  {
    parameters.get_interacting_bands().resize(0);
    
    parameters.get_interacting_bands().push_back(0);    
    parameters.get_interacting_bands().push_back(1);
  }

  LDA_grid_size() = parameters.get_H_k_grid_size();

  if(int(LDA_grid_size().size()) != DIMENSION)
    throw std::logic_error(__FUNCTION__);
}

template<typename DCA_point_group_type>
std::vector<int>& bilayer_model<DCA_point_group_type>::LDA_grid_size()
{
  static std::vector<int> v(0);
  return v;
}

template<typename DCA_point_group_type>
double* bilayer_model<DCA_point_group_type>::get_r_DCA_basis()
{
  static double* r_DCA = initialize_r_DCA_basis();
  return r_DCA;
}

template<typename DCA_point_group_type>
double* bilayer_model<DCA_point_group_type>::get_k_DCA_basis()
{
  static double* k_DCA = initialize_k_DCA_basis();
  return k_DCA;
}

template<typename DCA_point_group_type>
double* bilayer_model<DCA_point_group_type>::get_r_LDA_basis()
{
  static double* r_LDA = initialize_r_LDA_basis();
  return r_LDA;
}

template<typename DCA_point_group_type>
double* bilayer_model<DCA_point_group_type>::get_k_LDA_basis()
{
  static double* k_LDA = initialize_k_LDA_basis();
  return k_LDA;
}


template<typename DCA_point_group_type>
double* bilayer_model<DCA_point_group_type>::initialize_r_DCA_basis()
{
  static double* r_DCA = new double[4];
  
  r_DCA[0] = 1.0;  r_DCA[1] = 0.0;
  r_DCA[2] = 0.0;  r_DCA[3] = 1.0;

  return r_DCA;
}

template<typename DCA_point_group_type>
double* bilayer_model<DCA_point_group_type>::initialize_k_DCA_basis()
{
  static double* k_DCA = new double[4];
  
  k_DCA[0] = 2*M_PI;  k_DCA[1] = 0.;
  k_DCA[2] = 0.;      k_DCA[3] = 2*M_PI;

  return k_DCA;
}

template<typename DCA_point_group_type>
double* bilayer_model<DCA_point_group_type>::initialize_r_LDA_basis()
{
  static double* r_LDA = new double[4];
  
  r_LDA[0] = 1.;  r_LDA[1] = 0.;
  r_LDA[2] = 0.;  r_LDA[3] = 1.;

  return r_LDA;
}

template<typename DCA_point_group_type>
double* bilayer_model<DCA_point_group_type>::initialize_k_LDA_basis()
{
  static double* k_LDA = new double[4];

  k_LDA[0] = 2.*M_PI;  k_LDA[1] = 0.;
  k_LDA[2] = 0.;       k_LDA[3] = 2.*M_PI;

  return k_LDA;
}

template<typename DCA_point_group_type>
std::vector<int> bilayer_model<DCA_point_group_type>::get_flavors()
{
  static std::vector<int> flavors(BANDS);

  flavors[0]=0;
  flavors[1]=1;

  return flavors;
}

template<typename DCA_point_group_type>
std::vector<std::vector<double> > bilayer_model<DCA_point_group_type>::get_a_vectors()
{
  static std::vector<std::vector<double> > a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));

  return a_vecs;
}

template<typename DCA_point_group_type>
std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > bilayer_model<DCA_point_group_type>::get_orbital_permutations()
{
  static std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > permutations(0);
  return permutations;
}

template<typename DCA_point_group_type>
template<class domain, class parameters_type>
void bilayer_model<DCA_point_group_type>::initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
								   parameters_type&            parameters)
{
  H_interaction = 0.;

  for(int b_ind1=0; b_ind1<BANDS; b_ind1++)
    for(int s_ind1=0; s_ind1<2; s_ind1++)
      for(int b_ind2=0; b_ind2<BANDS; b_ind2++)
	for(int s_ind2=0; s_ind2<2; s_ind2++)
	  {
	    if(b_ind1==b_ind2 and s_ind1!=s_ind2)
	      H_interaction(b_ind1, s_ind1, b_ind2, s_ind2, 0) = U;

	    if(b_ind1!=b_ind2 and s_ind1!=s_ind2)
	      H_interaction(b_ind1, s_ind1, b_ind2, s_ind2, 0) = V;
	  }
}

template<typename DCA_point_group_type>
template<class domain>
void bilayer_model<DCA_point_group_type>::initialize_H_symmetries(FUNC_LIB::function<int, domain>& H_symmetries)
{
  H_symmetries = -1;

  H_symmetries(0,0, 0,0) = 0;
  H_symmetries(0,1, 0,1) = 0;

  H_symmetries(1,0, 1,0) = 1;
  H_symmetries(1,1, 1,1) = 1;
}

template<typename DCA_point_group_type>
template<class domain, class parameters_type>
void bilayer_model<DCA_point_group_type>::initialize_H_LDA(FUNC_LIB::function<std::complex<double>, domain>& H_LDA,
							   parameters_type&                        parameters)
{
  typedef typename parameters_type::b b;
  typedef typename parameters_type::s s;

  typedef typename parameters_type::k_LDA k_LDA; 

  std::vector<double> k;
  
  for(int k_ind=0; k_ind<k_LDA::dmn_size(); k_ind++)
    {
      k = k_LDA::parameter_type::get_elements()[k_ind];
      
      for(int b_ind1=0; b_ind1<b::dmn_size(); b_ind1++)
	for(int s_ind1=0; s_ind1<s::dmn_size(); s_ind1++)
	  for(int b_ind2=0; b_ind2<b::dmn_size(); b_ind2++)
	    for(int s_ind2=0; s_ind2<s::dmn_size(); s_ind2++)
	      H_LDA(b_ind1, s_ind1, b_ind2, s_ind2, k_ind) = get_LDA_Hamiltonians(parameters, k, b_ind1, s_ind1, b_ind2, s_ind2);    
    }  
}

template<typename DCA_point_group_type>
template<class parameters_type>
std::complex<double> bilayer_model<DCA_point_group_type>::get_LDA_Hamiltonians(parameters_type& parameters,
											  std::vector<double> k, 
											  int b1, int s1, int b2, int s2)
{
  const static std::complex<double> I(0,1);
  
  std::complex<double> H_LDA = 0.;

  if(s1==s2)
    {
      if(b1==b2)
	H_LDA += -2.*(t_x*cos(k[0])+t_y*cos(k[1]));
      
      if(b1!=b2)
	H_LDA += -t_z;
    }

  return H_LDA;
}

#endif
