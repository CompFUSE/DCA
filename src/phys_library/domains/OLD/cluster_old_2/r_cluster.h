//-*-C++-*-

#ifndef R_CLUSTER_H
#define R_CLUSTER_H

/*!
 *  \author Peter Staar
 */
template<class base_cluster_type>
class r_cluster<FULL, base_cluster_type>
{

public:

  typedef          r_cluster<FULL, base_cluster_type > this_type;
  typedef          k_cluster<FULL, base_cluster_type > dual_type;

//   typedef          base_cluster_type                   base_cluster;

//   typedef typename base_cluster_type::parameter_type                  parameter_type;
//   typedef typename parameter_type::point_group                        point_group_type;

//   typedef dmn_0<electron_band_domain>                                   b_dmn_t;
//   typedef dmn_0<this_type>                                              r_dmn_t;

//   typedef dmn_0<point_group_symmetry_domain<UNIT_CELL, base_cluster> >  sym_unit_cell_dmn_t;
//   typedef dmn_0<point_group_symmetry_domain<SUPER_CELL, base_cluster> > sym_super_cell_dmn_t;

  typedef std::vector<double>                                         element_type;

  const static cluster_representation_type cluster_representation = FULL;

//   const static int DIMENSION = parameter_type::DIMENSION;
  const static int DIMENSION = base_cluster_type::DIMENSION;

  typedef MATH_ALGORITHMS::harmonic_dmn_nD_type dmn_specifications_type;

public:

  static int                                get_size(); 

  static std::vector<std::vector<double> >& get_basis_vectors();
  static std::vector<std::vector<double> >& get_super_basis_vectors();

  static std::vector<std::vector<double> >& get_elements();

  

  static int origin_index();
 
//   static void interpolation_vector_to(std::vector<double>&                q, 
// 				      std::vector<std::complex<double> >& T_matrix);

//   static double find_minimal_distance(std::vector<double> vec_0,
// 				      std::vector<double> vec_1);

//   static std::pair<bool, std::vector<std::vector<double> > > find_equivalent_vectors_with_minimal_distance_to_origin(std::vector<double> R_vec);


  static int add     (int r_i, int r_j);
  static int subtract(int r_i, int r_j);

//   static std::vector<double> back_inside_lattice(std::vector<double> r);
//   static std::vector<double> back_inside_cluster(std::vector<double> r);

//   static std::vector<double> get_cluster_affine_coordinate(std::vector<double>& r);
//   static std::vector<double> get_lattice_affine_coordinate(std::vector<double>& r);

//   static function<std::pair<int,int>, dmn_2< dmn_2<r_dmn_t,b_dmn_t>, sym_super_cell_dmn_t > >& get_symmetry_matrix(){
//     //     static function<std::pair<int,int>, dmn_2< dmn_2<r_dmn_t,b_dmn_t>, sym_super_cell_dmn_t> > symmetry_matrix("r_symmetry_matrix_super_cell");
//     //     return symmetry_matrix;
    
//     return cluster_symmetry<this_type>::get_symmetry_matrix();
//   }

  static void print();

  template<class stream_type>
  static void to_JSON(stream_type& ss);

private :

  /*
  static std::vector<std::vector<double> >& get_sorted_elements();

  static int get_R_index(std::vector<double> r);

  static std::vector<double> back_inside_cluster_1D(std::vector<double>& r);
  static std::vector<double> back_inside_cluster_2D(std::vector<double>& r);
  static std::vector<double> back_inside_cluster_3D(std::vector<double>& r);

  static std::vector<element_type> initialize_sorted_elements();
  */

  static std::vector<int> initialize_subtract_matrix();
  static std::vector<int> initialize_addition_matrix();

  //static bool test_back_inside_cluster();
  //static bool test_subtract_matrix();
};


template<class base_cluster_type >
int r_cluster<FULL, base_cluster_type >::get_size()
{
  switch(cluster_representation)
    {
//     case IRREDUCIBLE:
//       return base_cluster_type::get_irreducible_r_cluster().size();
    case FULL:
      return base_cluster_type::get_cluster_size();
    default:
      cout << "ERROR in " << __PRETTY_FUNCTION__ << endl;
      throw std::logic_error(__FUNCTION__);
    }

  return -1;
}

template<class base_cluster_type >
std::vector<std::vector<double> >& r_cluster<FULL, base_cluster_type >::get_basis_vectors()
{
  return base_cluster_type::get_r_basis();
}

template<class base_cluster_type >
std::vector<std::vector<double> >& r_cluster<FULL, base_cluster_type >::get_super_basis_vectors()
{
  return base_cluster_type::get_r_super_basis();
}

template<class base_cluster_type >
std::vector<std::vector<double> >& r_cluster<FULL, base_cluster_type >::get_elements()
{
  switch(FULL)
    {
//     case IRREDUCIBLE:
//       return base_cluster_type::get_irreducible_r_cluster();
    case FULL:
      return base_cluster_type::get_r_cluster();
    default:
      {
	throw std::logic_error(__FUNCTION__);
	return base_cluster_type::get_r_cluster();
      }
    }
}


// template<class base_cluster_type >
// std::vector<std::vector<double> >& r_cluster<FULL, base_cluster_type >::get_sorted_elements()
// {
//   static std::vector<std::vector<double> > sorted_elements = initialize_sorted_elements();
//   return sorted_elements;
// }

// template<class base_cluster_type >
// std::vector<std::vector<double> > r_cluster<FULL, base_cluster_type >::initialize_sorted_elements()
// {
//   std::vector<std::vector<double> > sorted_elements = get_elements();

//   for(size_t l=0; l<sorted_elements.size(); l++)
//     sorted_elements[l].push_back(l);

//   sort(sorted_elements.begin(), sorted_elements.end(), VECTOR_OPERATIONS::IS_LARGER_VECTOR);

//   return sorted_elements;
// }

template<class base_cluster_type >
int r_cluster<FULL, base_cluster_type >::origin_index()
{
//   int index = -1;

//   for(int l=0; l<get_size(); l++)
//     if(VECTOR_OPERATIONS::L2_NORM(base_cluster_type::get_r_cluster()[l]) < 1.e-6)
//       index = l;

//   assert(index>-1 && index<get_size());

  static int index = cluster_operations::origin_index(get_elements());
  assert(VECTOR_OPERATIONS::L2_NORM(get_elements()[index])<1.e-6);

  return index;
}

// template<class base_cluster_type>
// int r_cluster<FULL, base_cluster_type >::get_R_index(std::vector<double> r)
// {
//   r.push_back(-1);

//   std::vector<std::vector<double> >& v = get_sorted_elements();

//   int index = int(lower_bound(v.begin(), v.end(), r, VECTOR_OPERATIONS::IS_LARGER_VECTOR) - v.begin());

//   index = v[index][DIMENSION];
//   assert(index>-1 and index<get_size());

//   {
//     r.pop_back();
//     assert(VECTOR_OPERATIONS::L2_NORM(r, get_elements()[index])<1.e-6);
//   }

//   return index;
// }


/*
template<class base_cluster_type >
void r_cluster<FULL, base_cluster_type >::interpolation_vector_to(std::vector<double>&                q_vec, 
										 std::vector<std::complex<double> >& T_matrix)
{
  assert(int(q_vec.size())    == DIMENSION);
  assert(int(T_matrix.size()) == base_cluster_type::get_cluster_size());

  for(int r=0; r<base_cluster_type::get_cluster_size(); r++)
    T_matrix[r] = 0.;

  std::vector<double> r_vec;
  std::vector<double> k_vec;

  for(int r=0; r<base_cluster_type::get_cluster_size(); r++){
    for(int k=0; k<base_cluster_type::get_cluster_size(); k++){
      
      r_vec = base_cluster_type::get_r_cluster()[r];
      k_vec = base_cluster_type::get_k_cluster()[k];
      
      double dot_prod = 0;
      for(size_t l=0; l<k_vec.size(); l++)
	dot_prod += k_vec[l]*(r_vec[l]-q_vec[l]);
      
      T_matrix[r] += std::complex<double>(cos(dot_prod), sin(dot_prod)); 
    }
  }

  for(int r=0; r<base_cluster_type::get_cluster_size(); r++)
    T_matrix[r] /= double(base_cluster_type::get_cluster_size());
}
*/

/*
template<class base_cluster_type >
double r_cluster<FULL, base_cluster_type >::find_minimal_distance(std::vector<double> vec_0,
										 std::vector<double> vec_1)
{
//   vec_0 = back_inside_cluster(vec_0);
//   vec_1 = back_inside_cluster(vec_1);

  vec_0 = cluster_operations::translate_inside_cluster(vec_0, get_super_basis());
  vec_1 = cluster_operations::translate_inside_cluster(vec_1, get_super_basis());

  double MIN_DISTANCE = VECTOR_OPERATIONS::L2_NORM(vec_0, vec_1);
  
  switch(DIMENSION)
    {
    case 2:
      
      for(int l0=-2; l0<=2; l0++){
	for(int l1=-2; l1<=2; l1++){
	  
	  std::vector<double> vec = vec_0;
	  
	  for(int d=0; d<DIMENSION; d++)
	    vec[d] += (l0*get_super_basis()[0][d] + l1*get_super_basis()[1][d]); 
	  
	  double distance = VECTOR_OPERATIONS::L2_NORM(vec, vec_1);
	  
	  if(distance < MIN_DISTANCE)
	    MIN_DISTANCE = distance;
	}
      }
      break;
      
    case 3:
      
      for(int l0=-2; l0<=2; l0++){
	for(int l1=-2; l1<=2; l1++){
	  for(int l2=-2; l2<=2; l2++){
	    
	    std::vector<double> vec = vec_0;
	    
	    for(int d=0; d<DIMENSION; d++)
	      vec[d] += (l0*get_super_basis()[0][d] + l1*get_super_basis()[1][d] + l2*get_super_basis()[2][d]); 
	    
	    double distance = VECTOR_OPERATIONS::L2_NORM(vec, vec_1);
	    
	    if(distance < MIN_DISTANCE)
	      MIN_DISTANCE = distance;
	  }
	}
      }
      break;
      
    default:
      throw std::logic_error(__FUNCTION__);
    }

    return MIN_DISTANCE;
}

template<class base_cluster_type >
std::pair<bool, std::vector<std::vector<double> > > r_cluster<FULL, base_cluster_type >::find_equivalent_vectors_with_minimal_distance_to_origin(std::vector<double> R_vec)
{
  const static double EPS          = 1.e-3;
  const static double ONE_PLUS_EPS = 1.+EPS;
  const static double ONE_MIN_EPS  = 1.-EPS;

  //R_vec = back_inside_cluster(R_vec);

  R_vec = cluster_operations::translate_inside_cluster(R_vec, get_super_basis());

  std::vector<double> origin(DIMENSION,0.);
  double MIN_DISTANCE = find_minimal_distance(origin, R_vec);
  
  bool IS_MINIMAL = VECTOR_OPERATIONS::L2_NORM(R_vec)>(MIN_DISTANCE+1.e-6)? false : true;

  std::vector<std::vector<double> > r_min;
 
  switch(DIMENSION)
    {
    case 2:

      for(int l0=-2; l0<=2; l0++){
	for(int l1=-2; l1<=2; l1++){
	    	    
	  std::vector<double> vec = R_vec;

	  for(int d=0; d<DIMENSION; d++)
	    vec[d] += (l0*get_super_basis()[0][d] + l1*get_super_basis()[1][d]); 

	  double distance = sqrt(VECTOR_OPERATIONS::L2_NORM(vec));

	  if(distance > sqrt(MIN_DISTANCE)*ONE_MIN_EPS-EPS && distance < sqrt(MIN_DISTANCE)*ONE_PLUS_EPS+EPS)
	    r_min.push_back(vec);
	}
      }
	      
      break;
      
    case 3:

      for(int l0=-2; l0<=2; l0++){
	for(int l1=-2; l1<=2; l1++){
	  for(int l2=-2; l2<=2; l2++){
	    
	    std::vector<double> vec = R_vec;
	    
	    for(int d=0; d<DIMENSION; d++)
	      vec[d] += (l0*get_super_basis()[0][d] + l1*get_super_basis()[1][d] + l2*get_super_basis()[2][d]); 
	    
	    double distance = sqrt(VECTOR_OPERATIONS::L2_NORM(vec));
	    
	    if(distance > sqrt(MIN_DISTANCE)*ONE_MIN_EPS-EPS && distance < sqrt(MIN_DISTANCE)*ONE_PLUS_EPS+EPS)
	      r_min.push_back(vec);
	  }
	}
      }
      break;
      
    default:
      throw std::logic_error(__FUNCTION__);
    }

  sort(r_min.begin(), r_min.end(), &VECTOR_OPERATIONS::SAME_VECTOR);
  int vec_size = unique(r_min.begin(), r_min.end(), &VECTOR_OPERATIONS::SAME_VECTOR)-r_min.begin();

  r_min.resize(vec_size);

  return std::pair<bool,  std::vector<std::vector<double> > >(IS_MINIMAL, r_min);
}
*/

template<class base_cluster_type >
int r_cluster<FULL, base_cluster_type >::add(int r_i, int r_j)
{
//   if(IS_EQUAL_TYPE<PsiMag_symmetry_2D, typename base_cluster_type::point_group_type>::check){
//     throw std::logic_error(__FUNCTION__);
// //     return get_PsiMag_crystal().add(r_j,r_i);
//   }
//   else
//     {
  static int              cluster_size = base_cluster_type::get_cluster_size();
  static std::vector<int> diff_indices = initialize_addition_matrix();
  
  assert(r_i >= 0 && r_i < cluster_size && r_j >= 0 && r_j < cluster_size);  
  return diff_indices[r_i + cluster_size*r_j];
//     }
}

template<class base_cluster_type >
int r_cluster<FULL, base_cluster_type >::subtract(int r_i, int r_j)
{
//   if(IS_EQUAL_TYPE<PsiMag_symmetry_2D, typename base_cluster_type::point_group_type>::check){
//     throw std::logic_error(__FUNCTION__);
// //     return get_PsiMag_crystal().subtract(r_j,r_i);
//   }
//   else
//     {
  static int              cluster_size = base_cluster_type::get_cluster_size();
  static std::vector<int> diff_indices = initialize_subtract_matrix();
  
  assert(r_i >= 0 && r_i < cluster_size && r_j >= 0 && r_j < cluster_size);  
  return diff_indices[r_i + cluster_size*r_j];
//     }
}

template<class base_cluster_type >
std::vector<int> r_cluster<FULL, base_cluster_type >::initialize_subtract_matrix()
{
  //cout << __PRETTY_FUNCTION__ << endl;

  //assert(test_back_inside_cluster());

  int cluster_size = base_cluster_type::get_cluster_size();
  std::vector<int> diff_indices(cluster_size*cluster_size);

  for(int i=0; i<cluster_size; i++)
    {
      for(int j=0; j<cluster_size; j++)
	{
	  std::vector<double> r_j_minus_r_i = base_cluster_type::get_r_cluster()[j];
	      
	  for(int l=0; l<DIMENSION; l++)
	    r_j_minus_r_i[l] -= base_cluster_type::get_r_cluster()[i][l];
	  
	  //r_j_minus_r_i  = back_inside_cluster(r_j_minus_r_i);
	  r_j_minus_r_i = cluster_operations::translate_inside_cluster(r_j_minus_r_i, get_super_basis_vectors());

	  /*
	  int index = -1;
	  for(int l=0; l<cluster_size; l++)
	    if(VECTOR_OPERATIONS::L2_NORM(r_j_minus_r_i, base_cluster_type::get_r_cluster()[l]) < 1.e-6)
	      index = l;
	  */
	  
	  //int index = get_R_index(r_j_minus_r_i);
	  int index = cluster_operations::index(r_j_minus_r_i, get_elements());

	  assert(index >= 0 && index < cluster_size);
	  diff_indices[i + j*cluster_size] = index;
	}
    }

  return diff_indices;
}

template<class base_cluster_type >
std::vector<int> r_cluster<FULL, base_cluster_type >::initialize_addition_matrix()
{
  //cout << __PRETTY_FUNCTION__ << endl;
  
  //assert(test_back_inside_cluster());

  int cluster_size = base_cluster_type::get_cluster_size();
  std::vector<int> diff_indices(cluster_size*cluster_size);

  for(int i=0; i<cluster_size; i++)
    {
      for(int j=0; j<cluster_size; j++)
	{
	  std::vector<double> r_j_plus_r_i = base_cluster_type::get_r_cluster()[j];
	  
	  for(int l=0; l<DIMENSION; l++)
	    r_j_plus_r_i[l] += base_cluster_type::get_r_cluster()[i][l];
	  
	  //r_j_plus_r_i  = back_inside_cluster(r_j_plus_r_i);
	  r_j_plus_r_i  = cluster_operations::translate_inside_cluster(r_j_plus_r_i, get_super_basis_vectors());

	  /*
	  int index = -1;
	  for(int l=0; l<cluster_size; l++)
	    if(VECTOR_OPERATIONS::L2_NORM(r_j_plus_r_i, base_cluster_type::get_r_cluster()[l]) < 1.e-6)
	      index = l;
	  */

	  //int index = get_R_index(r_j_plus_r_i);
	  int index = cluster_operations::index(r_j_plus_r_i, get_elements());

	  assert(index >= 0 && index < cluster_size);
	  diff_indices[i + j*cluster_size] = index;
	}
    }

  return diff_indices;
}


/*
template<class base_cluster_type >
std::vector<double> r_cluster<FULL, base_cluster_type >::back_inside_lattice(std::vector<double> r)
{

//   const static std::vector<std::vector<double> >& r_basis = base_cluster_type::get_r_basis();

//   std::vector<double> r_affine = get_lattice_affine_coordinate(r);

//   for(int d=0; d<DIMENSION; ++d){
//     while(r_affine[d]<-1.e-6)
//       r_affine[d] += 1.;
    
//     while(r_affine[d]>1-1.e-6)
//       r_affine[d] -= 1.;
//   }

//   std::vector<double> r_vec(DIMENSION,0);

//   for(int d1=0; d1<DIMENSION; ++d1)
//     for(int d0=0; d0<DIMENSION; ++d0)
//       r_vec[d0] += r_basis[d1][d0]*r_affine[d1];

//   return r_vec;


  std::vector<double> r_vec = cluster_operations::translate_inside_cluster(r, get_basis());

  return r_vec;
}
*/
 /*
template<class base_cluster_type >
std::vector<double> r_cluster<FULL, base_cluster_type >::back_inside_cluster(std::vector<double> r)
{
//   const static std::vector<std::vector<double> >& r_basis = base_cluster_type::get_r_super_basis();

//   std::vector<double> r_affine = get_cluster_affine_coordinate(r);

//   for(int d=0; d<DIMENSION; ++d){
//     while(r_affine[d]<-1.e-6)
//       r_affine[d] += 1.;
    
//     while(r_affine[d]>1-1.e-6)
//       r_affine[d] -= 1.;
//   }

//   std::vector<double> r_vec(DIMENSION,0);

//   for(int d1=0; d1<DIMENSION; ++d1)
//     for(int d0=0; d0<DIMENSION; ++d0)
//       r_vec[d0] += r_basis[d1][d0]*r_affine[d1];

//       return r_vec;

  std::vector<double> r_vec = cluster_operations::translate_inside_cluster(r, get_super_basis());

  return r_vec;
}
*/

  /* 
template<class base_cluster_type >
inline std::vector<double> r_cluster<FULL, base_cluster_type >::get_cluster_affine_coordinate(std::vector<double>& r)
{
//   static std::vector<std::vector<double> >& k_basis = base_cluster_type::get_k_super_basis();

//   static double one_over_two_pi = 1./(2.*M_PI);

//   std::vector<double> r_affine(DIMENSION,0);
  
//   for(int l1=0; l1<DIMENSION; l1++)
//     for(int l2=0; l2<DIMENSION; l2++)
//       r_affine[l1] += k_basis[l1][l2]*r[l2];

//   for(int l1=0; l1<DIMENSION; l1++)
//     r_affine[l1] *= one_over_two_pi;

//   return r_affine;

  std::vector<double> r_affine = VECTOR_OPERATIONS::COORDINATES(r, get_super_basis());
  return r_affine;
}
  */

   /*
template<class base_cluster_type >
inline std::vector<double> r_cluster<FULL, base_cluster_type >::get_lattice_affine_coordinate(std::vector<double>& r)
{
//   static std::vector<std::vector<double> >& k_basis = base_cluster_type::get_k_basis();

//   static double one_over_two_pi = 1./(2.*M_PI);

//   std::vector<double> r_affine(DIMENSION,0);
  
//   for(int l1=0; l1<DIMENSION; l1++)
//     for(int l2=0; l2<DIMENSION; l2++)
//       r_affine[l1] += k_basis[l1][l2]*r[l2];

//   for(int l1=0; l1<DIMENSION; l1++)
//     r_affine[l1] *= one_over_two_pi;

//   return r_affine;

  std::vector<double> r_affine = VECTOR_OPERATIONS::COORDINATES(r, get_basis());
  return r_affine;
}
   */





/*
template<class base_cluster_type >
bool r_cluster<FULL, base_cluster_type >::test_back_inside_cluster()
{
  static bool passed_test = false;

  if(!passed_test)
    {
      for(size_t l=0; l<base_cluster_type::get_r_cluster().size(); l++)
	{
	  std::vector<double> r1 = base_cluster_type::get_r_cluster()[l];
	  std::vector<double> r2 = back_inside_cluster(r1);

	  //cout << l << "\t" << VECTOR_OPERATIONS::L2_NORM(r1,r2) << "\n";
	  if(VECTOR_OPERATIONS::L2_NORM(r1,r2) > 1.e-6)
	    throw std::logic_error(__FUNCTION__);
	}
	
      passed_test = true;
    }

  return passed_test;
}
*/

// template<class base_cluster_type >
// bool r_cluster<FULL, base_cluster_type >::test_subtract_matrix()
// {
//   return true;
//   static bool passed_test = false;

//   if(!passed_test)
//     {
//       cout << __FUNCTION__ << endl;

//       std::vector<int> diff_indices = initialize_subtract_matrix();
//       int cluster_size = base_cluster_type::get_r_cluster().size();

//       for(size_t i=0; i<base_cluster_type::get_r_cluster().size(); i++)
// 	{
// 	  for(size_t j=0; j<base_cluster_type::get_r_cluster().size(); j++)
// 	    {
// // 	      if(diff_indices[i + cluster_size*j] != get_PsiMag_crystal().subtract(j,i))
// // 		throw std::logic_error(__FUNCTION__);
// 	    }
// 	}
	
//       passed_test = true;
//     }

//   return passed_test;
// }


template<class base_cluster_type >
template<class stream_type>
void r_cluster<FULL, base_cluster_type >::to_JSON(stream_type& ss)
{
  ss << "\"R-CLUSTER-OPERATIONS\" : \n {";

  ss << "\"subtract_ij = r_j-r_i\" : [\n";
  for(int i=0; i<get_size(); i++){
    ss << "[";
    for(int j=0; j<get_size(); j++){
      ss << subtract(i,j);
      j == get_size()-1 ? ss << "]": ss << ", ";
    }
    i == get_size()-1 ? ss << "\n" : ss << ",\n"; 
  }
  ss << "],\n"; 

  ss << "\"add_ij = r_j+r_i\" : [\n";
  for(int i=0; i<get_size(); i++){
    ss << "[";
    for(int j=0; j<get_size(); j++){
      ss << add(i,j);
      j == get_size()-1 ? ss << "]": ss << ", ";
    }
    i == get_size()-1 ? ss << "\n" : ss << ",\n"; 
  }

  ss << "]\n";

  ss << "}";
}

#endif
