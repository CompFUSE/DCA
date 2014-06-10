//-*-C++-*-

#ifndef SEARCH_MAXIMAL_SYMMETRY_GROUP_H
#define SEARCH_MAXIMAL_SYMMETRY_GROUP_H

/*!
 *  \ingroup SYMMETRIES
 *
 *  \author  Peter Staar
 *
 *  The operation is lattice compatible if there exists a integer matrix I, such that
 *
 *  O*V = V*I
 *
 *  The operation is lattice compatible if there exists a integer matrix I and a permutation matrix P, such that
 *
 *  O*a = a*P
 *
 */
template<class base_cluster_type, class point_group, symmetry_group_level_type symmetry_group_level, int INDEX>
class search_symmetry_group
{
private:

  typedef typename TypeAt<point_group, INDEX>::Result symmetry_type;
  typedef typename symmetry_type::base_type           group_action_type; 

  const static int DIMENSION = base_cluster_type::DIMENSION;

  typedef electron_band_domain               b_dmn_t;

//   typedef r_cluster<FULL, base_cluster_type> r_cluster_type;
//   typedef k_cluster<FULL, base_cluster_type> k_cluster_type;

  typedef typename base_cluster_type::r_cluster_type r_cluster_type;
  typedef typename base_cluster_type::k_cluster_type k_cluster_type;

  typedef point_group_symmetry_domain<symmetry_group_level, base_cluster_type> sym_dmn_t;

  typedef typename sym_dmn_t::element_type symmetry_element_type;

private:

  template<symmetry_group_level_type my_symmetry_group_level>
  struct intialize_T_matrices
  {
    static void execute(double* T, double* T_inv)
    {
      switch(my_symmetry_group_level)
	{
	case UNIT_CELL:
	  for(int i=0; i<DIMENSION; ++i)
	    for(int j=0; j<DIMENSION; ++j)
	      T[j+DIMENSION*i] = r_cluster_type::get_basis_vectors()[i][j];
	  break;
	  
	case SUPER_CELL:
	  for(int i=0; i<DIMENSION; ++i)
	    for(int j=0; j<DIMENSION; ++j)
	      T[j+DIMENSION*i] = r_cluster_type::get_super_basis_vectors()[i][j];
	  break;
	  
	default:
	  throw std::logic_error(__FUNCTION__);
	}
      
      {
	invert_plan<double> invert_pln(DIMENSION);
	
	memcpy(invert_pln.Matrix, &T[0]                     , sizeof(double)*square(DIMENSION));
	invert_pln.execute_plan();
	memcpy(&T_inv[0]        , invert_pln.inverted_matrix, sizeof(double)*square(DIMENSION));
      }
    }
  };

  static void back_inside_cluster(double* T_inv, double* v, double* w)
  {
    for(int i=0; i<DIMENSION; ++i)
      w[i] = 0.;

    for(int i=0; i<DIMENSION; i++)
      for(int j=0; j<DIMENSION; j++)
	w[i] += T_inv[i+DIMENSION*j]*v[j];
    
    for(int i=0; i<DIMENSION; ++i)
      {
	while(w[i]<-1.e-6)
	  w[i] += 1.;

	while(w[i]>1-1.e-6)
	  w[i] -= 1.;
      }
  }

  static bool is_duplicate()
  {
    if(sym_dmn_t::get_size() == 0)
      return false;
    else
      {
	bool is_a_duplicate = false;

	for(int l=0; l<sym_dmn_t::get_size(); l++){

	  bool is_this_a_duplicate = true;

	  for(int i=0; i<DIMENSION; i++)
	    for(int j=0; j<DIMENSION; j++)
	      if(fabs(sym_dmn_t::get_elements()[l].O[i+j*DIMENSION] - symmetry_type::matrix()[i+j*DIMENSION]) > 1.e-6)
		is_this_a_duplicate = false;

	  if(is_this_a_duplicate)
	    is_a_duplicate = true;
	}

	return is_a_duplicate;
      }
  }

  static bool is_lattice_compatible(double* T, double* T_inv)
  {
    double* permutation = new double[DIMENSION*DIMENSION];
    
    for(int i=0; i<DIMENSION; i++)
      for(int j=0; j<DIMENSION; j++)
	permutation[i+DIMENSION*j] = 0.;

    for(int i=0; i<DIMENSION; i++)
      for(int j=0; j<DIMENSION; j++)
	for(int l=0; l<DIMENSION; l++)
	  permutation[i+DIMENSION*j] += symmetry_type::matrix()[i+DIMENSION*l]*T[l+DIMENSION*j];

    bool is_compatible = true;

    std::vector<double> v  (DIMENSION, 0);
    std::vector<double> v_c(DIMENSION, 0);

    for(int i=0; i<DIMENSION; i++){

      for(int j=0; j<DIMENSION; j++)
	v[j] = permutation[j+DIMENSION*i];
      
      back_inside_cluster(&T_inv[0], &v[0], &v_c[0]);

      if(VECTOR_OPERATIONS::L2_NORM(v_c) > 1.e-6)
	is_compatible = false;
    }

    delete [] permutation;

    return is_compatible;
  }

  static bool is_unit_cell_compatible(double* T, double* T_inv)
  {
    int Na = b_dmn_t::get_size();

    double* permutation = new double[Na*DIMENSION];

    for(int j=0; j<Na; j++)
      for(int i=0; i<DIMENSION; i++)
	permutation[i+DIMENSION*j] = 0;//b_dmn_t::get_elements()[i].a_vec[j];

    for(int j=0; j<Na; j++)
      for(int i=0; i<DIMENSION; i++)
	for(int l=0; l<DIMENSION; l++)
	  permutation[i+DIMENSION*j] += symmetry_type::matrix()[i+DIMENSION*l]*b_dmn_t::get_elements()[j].a_vec[l];//T[l+DIMENSION*j];
    
    bool is_compatible = false;
    
    std::vector<double> v  (DIMENSION, 0);
    std::vector<double> v_c(DIMENSION, 0);
    
    for(int l=0; l<Na; l++){
      
      for(int j=0; j<DIMENSION; j++)
	v[j] = permutation[j+DIMENSION*l];
      
      // v_c is affine coordinates
      back_inside_cluster(&T_inv[0], &v[0], &v_c[0]);

      for(int j=0; j<DIMENSION; j++)
	v[j] = 0;//permutation[j+DIMENSION*i];

      for(int i=0; i<DIMENSION; i++)
	for(int j=0; j<DIMENSION; j++)
	  v[i] += symmetry_type::matrix()[i+DIMENSION*j]*v_c[j];

      for(int j=0; j<Na; j++)
	if(VECTOR_OPERATIONS::L2_NORM(v, b_dmn_t::get_elements()[j].a_vec) < 1.e-6
	   and b_dmn_t::get_elements()[l].flavor == b_dmn_t::get_elements()[j].flavor)
	  is_compatible = true;
    }

    delete [] permutation;

    return is_compatible;
  }

public:

  static void execute()
  {
    double* T     = new double[DIMENSION*DIMENSION];
    double* T_inv = new double[DIMENSION*DIMENSION];

    double* T_cell     = new double[DIMENSION*DIMENSION];
    double* T_cell_inv = new double[DIMENSION*DIMENSION];

    intialize_T_matrices<UNIT_CELL           >::execute(T_cell, T_cell_inv);
    intialize_T_matrices<symmetry_group_level>::execute(T     , T_inv);

//     {
//       cout << INDEX << endl;
//       for(int i=0; i<DIMENSION; i++){
// 	for(int j=0; j<DIMENSION; j++)
// 	  cout << symmetry_type::matrix()[i+j*DIMENSION] << "\t";
// 	cout << endl;
//       }
//       cout << endl;
//     }

    bool is_a_duplicate     = is_duplicate();

    bool unit_cell_compatible = is_unit_cell_compatible(T_cell, T_cell_inv);
    bool lattice_compatible   = is_lattice_compatible  (T     , T_inv);

    if(lattice_compatible && unit_cell_compatible && not is_a_duplicate){

      symmetry_element_type symmetry_element(DIMENSION);

      for(int i=0; i<DIMENSION; i++)
	for(int j=0; j<DIMENSION; j++)
	  symmetry_element.O[i+j*DIMENSION] = symmetry_type::matrix()[i+j*DIMENSION];

      for(int i=0; i<DIMENSION; i++)
	symmetry_element.t[i] = 0.;

      sym_dmn_t::get_elements().push_back(symmetry_element);

      sym_dmn_t::get_size() += 1;
    }

    delete [] T;
    delete [] T_inv;
    
    delete [] T_cell;
    delete [] T_cell_inv;

    search_symmetry_group<base_cluster_type, point_group, symmetry_group_level, INDEX-1>::execute();
  }
};

template<class cluster_type, class point_group, symmetry_group_level_type symmetry_group_level>
class search_symmetry_group<cluster_type, point_group, symmetry_group_level, -1>
{
public:
  static void execute()
  {}
};



/*!
 *  \ingroup SYMMETRIES
 *
 *  \author  Peter Staar
 */

template<class base_cluster_type, class point_group, symmetry_group_level_type symmetry_group_level>
class search_maximal_symmetry_group
{
public:

  typedef typename point_group::point_group_type_list point_group_type_list;

  const static int DIMENSION = base_cluster_type::DIMENSION;
  const static int MAX_SIZE  = Length<point_group_type_list>::value;

//   typedef r_cluster<FULL, base_cluster_type> r_cluster_type;
//   typedef k_cluster<FULL, base_cluster_type> k_cluster_type;

  typedef point_group_symmetry_domain<symmetry_group_level, base_cluster_type> sym_dmn_t;

public:

  static void execute()
  {
    sym_dmn_t::DIMENSION  = DIMENSION;
    sym_dmn_t::get_size() = 0;
    
    search_symmetry_group<base_cluster_type, point_group_type_list, symmetry_group_level, MAX_SIZE-1>::execute();
   
//     sym_dmn_t::to_JSON(std::cout);
//     throw std::logic_error(__FUNCTION__);
  }

};


#endif
