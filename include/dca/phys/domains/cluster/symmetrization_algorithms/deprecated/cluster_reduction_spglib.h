//-*-C++-*-

#ifndef CLUSTER_REDUCTION_SPG_LIB_H
#define CLUSTER_REDUCTION_SPG_LIB_H

namespace SPG_LIB
{
  extern "C"
  {
#include "spglib.h"
  }
}

#include "search_maximal_symmetry_group_spglib.h"

/*!
 *  \author Peter Staar
 *
 *
 *   r_{cart} = T \times r_{aff}
 *
 *   r'_{aff} = O \times r_{aff} + t
 *   r'_{cart} = T O T^{-1} r'_{cart} + T t
 */
template<class base_cluster_type>
class cluster_reduction<base_cluster_type, symmetry_package_spg_lib>
{
public:

  const static int DIMENSION = base_cluster_type::DIMENSION;

  const static int MAX_R_a = 1024;
  const static int MAX_OPS = 1024;

  typedef r_cluster<FULL, base_cluster_type> r_cluster_type;
  typedef k_cluster<FULL, base_cluster_type> k_cluster_type;

  typedef dmn_0<r_cluster_type> r_dmn_t;
  typedef dmn_0<k_cluster_type> k_dmn_t;

  typedef dmn_0<electron_band_domain>        b_dmn_t;

  typedef point_group_symmetry_domain<UNIT_CELL , base_cluster_type> sym_unit_cell_t;
  typedef point_group_symmetry_domain<SUPER_CELL, base_cluster_type> sym_super_cell_t;

  typedef dmn_0<sym_unit_cell_t>  sym_unit_cell_dmn_t;
  typedef dmn_0<sym_super_cell_t> sym_super_cell_dmn_t;
  
public:
  
  cluster_reduction();
  ~cluster_reduction();
  
  void execute();

private:

  void generate_transformation_matrices();

  void generate_cartesian_cluster();
  void generate_affine_cluster();

  void generate_permutations();
  void reduce_unit_cell_symmetries_to_super_cell_symmetries();

  void set_symmetry_matrix();

  void set_r_symmetry_matrix();
  void set_k_symmetry_matrix();

  int find_k_index(std::vector<double> k);

  void print_on_shell();

private:
  
  double b[3][3];
  double B[3][3];

  int    flavors [MAX_R_a];

  double a_i_cart[MAX_R_a][3];
  double a_i_frac[MAX_R_a][3];

  double r_i_cart[MAX_R_a][3];

  std::vector<std::vector<int> > permutations;
};

template<class base_cluster_type>
cluster_reduction<base_cluster_type, symmetry_package_spg_lib>::cluster_reduction():
  permutations(0)
{}
 
template<class base_cluster_type>
cluster_reduction<base_cluster_type, symmetry_package_spg_lib>::~cluster_reduction()
{}

template<class base_cluster_type>
void cluster_reduction<base_cluster_type, symmetry_package_spg_lib>::execute()
{
  generate_transformation_matrices();

  generate_cartesian_cluster();

  generate_permutations();

  {
    int Na = b_dmn_t::dmn_size();

    search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, UNIT_CELL> search_maximal_symmetry_group(DIMENSION, b, b, 1, r_i_cart, Na, flavors, a_i_cart);
    search_maximal_symmetry_group.initialize_symmetry_domain();

//     sym_unit_cell_t::to_JSON(std::cout);
  }

  {
    int Na = b_dmn_t::dmn_size();
    int Ns = r_dmn_t::dmn_size();

    search_maximal_symmetry_group<base_cluster_type, symmetry_package_spg_lib, SUPER_CELL> search_maximal_symmetry_group(DIMENSION, b, B, Ns, r_i_cart, Na, flavors, a_i_cart);
    search_maximal_symmetry_group.initialize_symmetry_domain();

//     sym_super_cell_t::to_JSON(std::cout);
  }

//   {
//     generate_permutations();
//     reduce_unit_cell_symmetries_to_super_cell_symmetries();
//   }

  set_symmetry_matrix();
}


template<class base_cluster_type>
void cluster_reduction<base_cluster_type, symmetry_package_spg_lib>::generate_transformation_matrices()
{
  for(int i=0; i<3; ++i){
    for(int j=0; j<3; ++j){
      b[i][j] = (i==j)? 1.:0.;
      B[i][j] = (i==j)? 1.:0.; 
    }
  }

  for(int i=0; i<DIMENSION; ++i){
    for(int j=0; j<DIMENSION; ++j){
      b[i][j] = r_cluster_type::get_basis()      [i][j];
      B[i][j] = r_cluster_type::get_super_basis()[i][j];
     }
   }
}

template<class base_cluster_type>
void cluster_reduction<base_cluster_type, symmetry_package_spg_lib>::generate_cartesian_cluster()
{
  int N_cell = r_dmn_t::dmn_size()*b_dmn_t::dmn_size();

  if(N_cell>MAX_R_a)
    throw std::logic_error(__FUNCTION__);

  {// unit-cell
    for(int j=0; j<b_dmn_t::dmn_size(); ++j){
      
      std::vector<double> t = b_dmn_t::get_elements()[j].a_vec;
      
      if(DIMENSION==2)
	t.push_back(0.);
      
      for(int l=0; l<3; ++l)
	a_i_cart[j][l] = t[l];
      
      flavors[j] = b_dmn_t::get_elements()[j].flavor;
    }
  }

  {// cluster
    for(int j=0; j<r_dmn_t::dmn_size(); ++j){
      
      std::vector<double> t = r_dmn_t::get_elements()[j];
      
      if(DIMENSION==2)
	t.push_back(0.);
      
      for(int l=0; l<3; ++l)
	r_i_cart[j][l] = t[l];
    }
  }
}

template<class base_cluster_type>
void cluster_reduction<base_cluster_type, symmetry_package_spg_lib>::generate_permutations()
{
  permutations.resize(0);
  
  std::vector<int> p(b_dmn_t::dmn_size(), 0);

  for(int b_i=0; b_i<b_dmn_t::dmn_size(); ++b_i)
    p[b_i] = b_i;

  do 
    {      
      bool is_good_permutation = true;
      
      for(int b_i=0; b_i<b_dmn_t::dmn_size(); ++b_i)
	{
	  std::vector<double> ai = b_dmn_t::get_elements()[  b_i ].a_vec;
	  std::vector<double> aj = b_dmn_t::get_elements()[p[b_i]].a_vec;
	  
	  double diff = VECTOR_OPERATIONS::L2_NORM(ai, aj);
	  
	  if(flavors[b_i] != flavors[p[b_i]] or diff>1.e-6)
	    is_good_permutation = false;
	}
      
      if(is_good_permutation)
	permutations.push_back(p);

//       for(int b_i=0; b_i<b_dmn_t::dmn_size(); ++b_i)
// 	cout << "\t" << p[b_i];
	
//       if(is_good_permutation)
// 	cout << "\tOK\n";
//       else
// 	cout << "\tNOT OK\n";
    }
  while (std::next_permutation(p.begin(),p.end()));
}

template<class base_cluster_type>
void cluster_reduction<base_cluster_type, symmetry_package_spg_lib>::reduce_unit_cell_symmetries_to_super_cell_symmetries()
{
  sym_super_cell_t::DIMENSION  = DIMENSION;
  sym_super_cell_t::get_size() = 0;

  std::vector<double> r_vec_i(DIMENSION, 0.);
  std::vector<double> r_vec_j(DIMENSION, 0.);

  std::vector<double> a_vec_i(DIMENSION, 0.);
  std::vector<double> a_vec_j(DIMENSION, 0.);

  std::vector<double> ri_plus_ai      (DIMENSION, 0.);
  std::vector<double> trafo_ri_plus_ai(DIMENSION, 0.);
  
  std::vector<double> rj_plus_aj    (DIMENSION, 0.);

  for(int l=0; l<sym_unit_cell_dmn_t::dmn_size(); ++l){
    for(size_t p=0; p<permutations.size(); ++p){

      bool is_super_cell_symmetry = true;
      
      for(int r_i=0; r_i<r_dmn_t::dmn_size(); ++r_i){
	for(int b_i=0; b_i<b_dmn_t::dmn_size(); ++b_i){
	  
	  r_vec_i = r_dmn_t::get_elements()[r_i];
	  a_vec_i = b_dmn_t::get_elements()[permutations[p][b_i]].a_vec;
	  
	  for(int d=0; d<DIMENSION; ++d)
	    ri_plus_ai[d] = r_vec_i[d] + a_vec_i[d];
	  
	  sym_unit_cell_dmn_t::get_elements()[l].transform(&ri_plus_ai[0], &trafo_ri_plus_ai[0]);	
	  trafo_ri_plus_ai = r_cluster_type::back_inside_cluster(trafo_ri_plus_ai);
	  
	  bool is_inside_cluster = false;
	  for(int r_j=0; r_j<r_dmn_t::dmn_size(); ++r_j){
	    for(int b_j=0; b_j<b_dmn_t::dmn_size(); ++b_j){
	      
	      r_vec_j = r_dmn_t::get_elements()[r_j];
	      a_vec_j = b_dmn_t::get_elements()[b_j].a_vec;
	      
	      for(int d=0; d<DIMENSION; ++d)
		rj_plus_aj[d] = r_vec_j[d] + a_vec_j[d];
	      
	      rj_plus_aj = r_cluster_type::back_inside_cluster(rj_plus_aj);
	      
	      if(VECTOR_OPERATIONS::L2_NORM(rj_plus_aj, trafo_ri_plus_ai)<1.e-6 and
		 b_dmn_t::get_elements()[permutations[p][b_i]].flavor == b_dmn_t::get_elements()[b_j].flavor)
		is_inside_cluster = true;
	    }
	  }
	  
	  if(not is_inside_cluster)
	    is_super_cell_symmetry = false;
	}
      }

      if(is_super_cell_symmetry){

	sym_super_cell_t::get_elements().push_back(sym_unit_cell_dmn_t::get_elements()[l]);

	sym_super_cell_t::get_elements().back().set_permutation(permutations[p]);

	sym_super_cell_t::get_size() += 1;
      }
    }
  }
}

template<class base_cluster_type>
void cluster_reduction<base_cluster_type, symmetry_package_spg_lib>::set_symmetry_matrix()
{
  set_r_symmetry_matrix();
  
  set_k_symmetry_matrix();
}

template<class base_cluster_type>
void cluster_reduction<base_cluster_type, symmetry_package_spg_lib>::set_r_symmetry_matrix()
{
  FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<r_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& symmetry_matrix = r_cluster_type::get_symmetry_matrix();

  for(int r_i=0; r_i<r_dmn_t::dmn_size(); ++r_i){
    for(int b_i=0; b_i<b_dmn_t::dmn_size(); ++b_i){

      for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l){
	

	//std::vector<int> permutation = sym_super_cell_dmn_t::get_elements()[l].get_permutation();

	symmetry_matrix(r_i,b_i, l) = std::pair<int,int>(-1, -1);

	std::vector<double> r_plus_a = VECTOR_OPERATIONS::ADD(r_dmn_t::get_elements()[r_i], 
							      b_dmn_t::get_elements()[b_i].a_vec);
	std::vector<double> trafo_r_plus_a(DIMENSION, 0);
	  
	sym_super_cell_dmn_t::get_elements()[l].transform(&r_plus_a[0], &trafo_r_plus_a[0]);
	
	trafo_r_plus_a = r_cluster_type::back_inside_cluster(trafo_r_plus_a);
	  
	for(int r_j=0; r_j<r_dmn_t::dmn_size(); ++r_j){
	  for(int b_j=0; b_j<b_dmn_t::dmn_size(); ++b_j){
	    
	    std::vector<double> rj_plus_aj = VECTOR_OPERATIONS::ADD(r_dmn_t::get_elements()[r_j], 
								    b_dmn_t::get_elements()[b_j].a_vec);

	    rj_plus_aj = r_cluster_type::back_inside_cluster(rj_plus_aj);
	      
	    if(VECTOR_OPERATIONS::L2_NORM(rj_plus_aj, trafo_r_plus_a)<1.e-6 and
	       b_dmn_t::get_elements()[b_i].flavor == b_dmn_t::get_elements()[b_j].flavor)
	      symmetry_matrix(r_i,b_i, l) = std::pair<int,int>(r_j,b_j);
	  }
	}
      }
    }
  }
}

template<class base_cluster_type>
void cluster_reduction<base_cluster_type, symmetry_package_spg_lib>::set_k_symmetry_matrix()
{
  FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<r_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& r_symmetry_matrix = r_cluster_type::get_symmetry_matrix();
  FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<k_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& k_symmetry_matrix = k_cluster_type::get_symmetry_matrix();

  for(int i=0; i<k_dmn_t::dmn_size(); ++i){
    for(int j=0; j<b_dmn_t::dmn_size(); ++j){

      for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l){

	std::vector<double> k       = k_dmn_t::get_elements()[i];
	std::vector<double> trafo_k(DIMENSION,0);

	sym_super_cell_dmn_t::get_elements()[l].linear_transform(&k[0], &trafo_k[0]);

	k_symmetry_matrix(i,j, l).first  = find_k_index(trafo_k);
	k_symmetry_matrix(i,j, l).second = r_symmetry_matrix(i, j, l).second;
      }
    }
  }
}

template<class cluster_type>
int cluster_reduction<cluster_type, symmetry_package_spg_lib>::find_k_index(std::vector<double> k)
{
  assert(k_cluster_type::test_back_inside_cluster());

  int index=-1;

  k = k_cluster_type::back_inside_cluster(k);

  for(int k_ind=0; k_ind<k_dmn_t::dmn_size(); ++k_ind)
    if(VECTOR_OPERATIONS::L2_NORM(k_dmn_t::get_elements()[k_ind], k)<1.e-6)
      index = k_ind;

  if(index<0)
    throw std::logic_error(__FUNCTION__);

  return index;
}


template<class cluster_type>
void cluster_reduction<cluster_type, symmetry_package_spg_lib>::print_on_shell()
{
  if(true)
    {
      FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<k_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& symmetry_matrix = k_cluster_type::get_symmetry_matrix();
      
      cout << "\n\t k-space symmetries : \n";
      for(int i=0; i<k_dmn_t::dmn_size(); ++i){
	for(int j=0; j<b_dmn_t::dmn_size(); ++j){
	  
	  cout << "\t" << i << ", " << j << "   (" << flavors[j] << ")\t|\t";
	  
	  for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l)
	    cout << "\t" << symmetry_matrix(i,j, l).first << ", " << symmetry_matrix(i,j, l).second;
	  
	  cout << "\n";
	}
      }
      cout << "\n";
    }

  if(true)
    {
      FUNC_LIB::function<std::pair<int,int>, dmn_2< dmn_2<r_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& symmetry_matrix = r_cluster_type::get_symmetry_matrix();
     
      cout << "\n\t r-space symmetries : \n";
      for(int i=0; i<r_dmn_t::dmn_size(); ++i){
	for(int j=0; j<b_dmn_t::dmn_size(); ++j){
	  
	  cout << "\t" << i << ", " << j << "   (" << flavors[j] << ")\t|\t";
	  
	  for(int l=0; l<sym_super_cell_dmn_t::dmn_size(); ++l)
	    cout << "\t" << symmetry_matrix(i,j, l).first << ", " << symmetry_matrix(i,j, l).second;
	  
	  cout << "\n";
	}
      }
      cout << "\n";
    }
}

#endif
