//-*-C++-*-

#ifndef TIGHT_BINDING_MODEL_H
#define TIGHT_BINDING_MODEL_H

/*!
 *   \author Peter Staar
 */
template<typename lattice_type, typename interaction_type>
class tight_binding_model
{
 public:

  typedef lattice_type     lattice;
  typedef interaction_type interaction;

public:

  // taken care off via parameters !
  static int&              get_DCA_size();
  static int&              get_LDA_size();

  static std::vector<int>& DCA_grid_size();
  static std::vector<int>& LDA_grid_size();

  // taken care of via lattice-type
  typedef typename lattice_type::LDA_point_group LDA_point_group;
  typedef typename lattice_type::DCA_point_group DCA_point_group;

  static const int DIMENSION = lattice_type::DIMENSION;
  static const int BANDS     = lattice_type::BANDS;

//   static const cluster_shape_type DCA_cluster_shape = lattice_type::DCA_cluster_shape;
//   static const cluster_shape_type LDA_cluster_shape = lattice_type::LDA_cluster_shape;

  static double* get_r_DCA_basis();
  static double* get_k_DCA_basis();

  static double* get_r_LDA_basis();
  static double* get_k_LDA_basis();

  static std::vector<int>                  get_flavors();
  static std::vector<std::vector<double> > get_a_vectors();

  static std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > get_orbital_permutations();

  template<class domain, class parameters_type>
  static void initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
				       parameters_type&            parameters);

  template<class domain>
  static void initialize_H_symmetries(FUNC_LIB::function<int, domain>& H_interactions_symmetries);

  template<class domain, class parameters_type>
  static void initialize_H_LDA(FUNC_LIB::function<std::complex<double> , domain >& H_LDA,
			       parameters_type&                          parameters);

  template<class parameters_type>
  static void initialize(parameters_type& parameters);

private:

  template<class parameters_type>
  static void initialize_Bett_cluster(parameters_type& parameters);

  template<class parameters_type>
  static void initialize_default(parameters_type& parameters);
};

template<typename lattice_type, typename interaction_type>
int& tight_binding_model<lattice_type, interaction_type>::get_DCA_size()
{
  static int DCA_size = -1;
  return DCA_size;
}

template<typename lattice_type, typename interaction_type>
int& tight_binding_model<lattice_type, interaction_type>::get_LDA_size()
{
  static int LDA_size = -1;
  return LDA_size;
}

template<typename lattice_type, typename interaction_type>
std::vector<int>& tight_binding_model<lattice_type, interaction_type>::DCA_grid_size()
{
  static std::vector<int> v(0);
  return v;
}

template<typename lattice_type, typename interaction_type>
std::vector<int>& tight_binding_model<lattice_type, interaction_type>::LDA_grid_size()
{
  static std::vector<int> v(0);
  return v;
}

template<typename lattice_type, typename interaction_type>
double* tight_binding_model<lattice_type, interaction_type>::get_r_DCA_basis()
{
  static double* r_DCA = lattice_type::initialize_r_DCA_basis();
  return r_DCA;
}

// template<typename lattice_type, typename interaction_type>
// double* tight_binding_model<lattice_type, interaction_type>::get_k_DCA_basis()
// {
//   static double* k_DCA = lattice_type::initialize_k_DCA_basis(); 
//   return k_DCA;
// }

template<typename lattice_type, typename interaction_type>
double* tight_binding_model<lattice_type, interaction_type>::get_r_LDA_basis()
{
  static double* r_LDA = lattice_type::initialize_r_LDA_basis();
  return r_LDA;
}

// template<typename lattice_type, typename interaction_type>
// double* tight_binding_model<lattice_type, interaction_type>::get_k_LDA_basis()
// {
//   static double* k_LDA = lattice_type::initialize_k_LDA_basis();
//   return k_LDA;
// }

template<typename lattice_type, typename interaction_type>
std::vector<int> tight_binding_model<lattice_type, interaction_type>::get_flavors()
{
  return lattice_type::get_flavors();
}

template<typename lattice_type, typename interaction_type>
std::vector<std::vector<double> > tight_binding_model<lattice_type, interaction_type>::get_a_vectors()
{
  return lattice_type::get_a_vectors();
}

/*
template<typename lattice_type, typename interaction_type>
std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > tight_binding_model<lattice_type, interaction_type>::get_orbital_permutations()
{
  static std::vector<std::pair<std::pair<int,int>, std::pair<int,int> > > permutations = lattice_type::get_orbital_permutations();
  return permutations;
}
*/

template<typename lattice_type, typename interaction_type>
template<class domain, class parameters_type>
void tight_binding_model<lattice_type, interaction_type>::initialize_H_interaction(FUNC_LIB::function<double , domain >& H_interaction,
										   parameters_type&            parameters)
{
  lattice_type::initialize_H_interaction(H_interaction, parameters);
}

template<typename lattice_type, typename interaction_type>
template<class domain>
void tight_binding_model<lattice_type, interaction_type>::initialize_H_symmetries(FUNC_LIB::function<int, domain>& H_symmetry)
{
  lattice_type::initialize_H_symmetry(H_symmetry);
}

template<typename lattice_type, typename interaction_type>
template<class domain, class parameters_type>
void tight_binding_model<lattice_type, interaction_type>::initialize_H_LDA(FUNC_LIB::function<std::complex<double> , domain >& H_LDA,
									   parameters_type&                          parameters)
{
  typedef typename parameters_type::k_LDA k_LDA; 
  typedef typename parameters_type::LDA_k_cluster_type LDA_k_cluster_type;
  typedef typename parameters_type::b b;
  typedef typename parameters_type::s s;

  std::vector<double> k;
  
  for(int k_ind=0; k_ind<k_LDA::dmn_size(); k_ind++)
    {
      k = LDA_k_cluster_type::get_elements()[k_ind];
      
      for(int b_ind1=0; b_ind1<b::dmn_size(); b_ind1++)
	for(int s_ind1=0; s_ind1<s::dmn_size(); s_ind1++)
	  for(int b_ind2=0; b_ind2<b::dmn_size(); b_ind2++)
	    for(int s_ind2=0; s_ind2<s::dmn_size(); s_ind2++)
	      H_LDA(b_ind1, s_ind1, b_ind2, s_ind2, k_ind) = lattice_type::get_LDA_Hamiltonians(parameters, k, b_ind1, s_ind1, b_ind2, s_ind2);    
    }  
}


template<typename lattice_type, typename interaction_type>
template<class parameters_type>
void tight_binding_model<lattice_type, interaction_type>::initialize(parameters_type& /*parameters*/)
{
/*
  // compile-time switch

    switch(DCA_cluster_shape)
      {
      case BETT_CLUSTER:
	initialize_Bett_cluster(parameters);
	break;
	
      case PARALLELEPIPED :
	initialize_default(parameters);
	break;
	
      default:
	throw std::logic_error(__FUNCTION__);
      }
*/
}

/*
template<typename lattice_type, typename interaction_type>
template<class parameters_type>
void tight_binding_model<lattice_type, interaction_type>::initialize_Bett_cluster(parameters_type& parameters)
{
  LDA_grid_size() = parameters.get_H_k_grid_size();//parameters.get_LDA_grid_size();

  if(int(LDA_grid_size().size()) != DIMENSION)
    throw std::logic_error(__FUNCTION__);

  get_LDA_size() = 1;
  for(int l=0; l<DIMENSION; l++)
    get_LDA_size() *= LDA_grid_size()[l];
}

template<typename lattice_type, typename interaction_type>
template<class parameters_type>
void tight_binding_model<lattice_type, interaction_type>::initialize_default(parameters_type& parameters)
{
  throw std::logic_error(__FUNCTION__);

//   DCA_grid_size() = parameters.get_DCA_grid_size();
//   LDA_grid_size() = parameters.get_LDA_grid_size();

//   if(int(DCA_grid_size().size()) != DIMENSION || int(LDA_grid_size().size()) != DIMENSION)
//     throw std::logic_error(__FUNCTION__);

//   get_DCA_size() = 1;
//   get_LDA_size() = 1;
//   for(int l=0; l<DIMENSION; l++){
//     get_DCA_size() *= DCA_grid_size()[l];
//     get_LDA_size() *= LDA_grid_size()[l];
//   }
}
*/

#endif
