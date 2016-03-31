//-*-C++-*-

#ifndef CLUSTER_OPERATIONS_H
#define CLUSTER_OPERATIONS_H

#include <utility>
#include <vector>
#include "cluster_typedefs.hpp"
#include "math_library/geometry_library/vector_operations/include_vector_operations.h"

/*!
 *  \author Peter Staar
 */
class cluster_operations
{
 public:
  
  template<typename scalar_type>
  static int index(std::vector<scalar_type>&               element,
		   std::vector<std::vector<scalar_type> >& elements,
		   CLUSTER_SHAPE                           shape);

  template<typename scalar_type>
  static int origin_index(std::vector<std::vector<scalar_type> >& elements,
			  CLUSTER_SHAPE                           shape);
  
  template<typename scalar_type>
  static std::vector<scalar_type> translate_inside_cluster(std::vector<scalar_type>& r,
							   std::vector<std::vector<scalar_type> >& basis);

  template<typename scalar_type>
  static bool test_translate_inside_cluster(std::vector<std::vector<scalar_type> >& elements,
					    std::vector<std::vector<scalar_type> >& basis);

  template<typename scalar_type>
  static scalar_type minimal_distance(std::vector<scalar_type> vec_0,
				      std::vector<scalar_type> vec_1,
				      std::vector<std::vector<scalar_type> >& basis);
  template<typename scalar_type>
  static bool is_minimal(std::vector<scalar_type> R_vec,
			 std::vector<std::vector<scalar_type> >& basis);

  template<typename scalar_type>
  static std::vector<std::vector<scalar_type> > equivalent_vectors(std::vector<scalar_type> R_vec,
								   std::vector<std::vector<scalar_type> >& basis);

  template<typename cluster_type, typename scalar_type>
  static std::pair<std::vector<scalar_type>, int> find_closest_cluster_vector(const std::vector<scalar_type>& input_vec, const double tol);
};


template<typename scalar_type>
int cluster_operations::index(std::vector<scalar_type>&               element,
			      std::vector<std::vector<scalar_type> >& sorted_elements,
			      CLUSTER_SHAPE                           shape)
{
  /*
  int index = lower_bound(sorted_elements.begin(), sorted_elements.end(), 
			  element, VECTOR_OPERATIONS::IS_LARGER_VECTOR<scalar_type>) - sorted_elements.begin();
  */

  int index = -1;

  switch(shape)
    {
    case BRILLOUIN_ZONE: // the k-vectors in the brillouin zone are sorted according to VECTOR_OPERATIONS::IS_LARGER_VECTOR !
      index = lower_bound(sorted_elements.begin(), sorted_elements.end(), 
			  element, VECTOR_OPERATIONS::IS_LARGER_VECTOR<scalar_type>) - sorted_elements.begin();
      break;

    case PARALLELLEPIPEDUM:
      index = find(sorted_elements.begin(), sorted_elements.end(), element) - sorted_elements.begin();
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }

  assert(index>-1 and index<sorted_elements.size());

  if(VECTOR_OPERATIONS::L2_NORM(element, sorted_elements[index])>1.e-6)
    {
      std::cout << "\n\t " << __FUNCTION__ << "\t" << index << "\n"; 
      VECTOR_OPERATIONS::PRINT(element); std::cout << "\n";  
      VECTOR_OPERATIONS::PRINT(sorted_elements[index]); std::cout << "\n";
      std::cout << "\n\n";
    }

  assert(VECTOR_OPERATIONS::L2_NORM(element, sorted_elements[index])<1.e-6);
  
  return index;
}

template<typename scalar_type>
int cluster_operations::origin_index(std::vector<std::vector<scalar_type> >& elements,
				     CLUSTER_SHAPE                           shape)
{
  std::vector<double> origin(elements[0].size(), 0.);
  
  return index(origin, elements, shape);
}

template<typename scalar_type>
std::vector<scalar_type> cluster_operations::translate_inside_cluster(std::vector<scalar_type>& r,
								      std::vector<std::vector<scalar_type> >& basis)
{
  int DIMENSION = r.size();
  
  std::vector<scalar_type> r_affine = VECTOR_OPERATIONS::COORDINATES(r, basis);
  
  for(size_t d=0; d<r.size(); d++)
    {
      while(r_affine[d]<-1.e-6)
	r_affine[d] += 1.;
      
      while(r_affine[d]>1-1.e-6)
	r_affine[d] -= 1.;
    }
  
  std::vector<scalar_type> r_vec(r.size(), 0.);
  
  for(size_t d1=0; d1<DIMENSION; ++d1)
    for(size_t d0=0; d0<DIMENSION; ++d0)
      r_vec[d0] += basis[d1][d0]*r_affine[d1];
  
  return r_vec;
}

template<typename scalar_type>
bool cluster_operations::test_translate_inside_cluster(std::vector<std::vector<scalar_type> >& elements,
						       std::vector<std::vector<scalar_type> >& basis)
{
  static bool passed_test = false;

  if(!passed_test)
    {
      std::vector<scalar_type> k1, k2;

      for(size_t l=0; l<elements.size(); l++)
	{
	  k1 = elements[l];
	  k2 = translate_inside_cluster(k1, basis);

	  if(VECTOR_OPERATIONS::L2_NORM(k1,k2) > 1.e-6)
	    throw std::logic_error(__FUNCTION__);
	}
	
      passed_test = true;
    }

  return passed_test;
}

template<typename scalar_type>
scalar_type cluster_operations::minimal_distance(std::vector<scalar_type> vec_0,
						 std::vector<scalar_type> vec_1,
						 std::vector<std::vector<scalar_type> >& basis)
{
  int DIMENSION = vec_0.size();

  vec_0 = translate_inside_cluster(vec_0, basis);
  vec_1 = translate_inside_cluster(vec_1, basis);

  scalar_type MIN_DISTANCE = VECTOR_OPERATIONS::L2_NORM(vec_0, vec_1);
  
  switch(DIMENSION)
    {
    case 1:
      {
	for(int l0=-1; l0<=1; l0++)
	  {
	    std::vector<scalar_type> vec = vec_0;
	  
	    for(int d=0; d<DIMENSION; d++)
	      vec[d] += (l0*basis[0][d]); 
	    
	    scalar_type distance = VECTOR_OPERATIONS::L2_NORM(vec, vec_1);
	    
	    if(distance < MIN_DISTANCE)
	      MIN_DISTANCE = distance;
	  }
      }
      break;

    case 2:
      {
	for(int l0=-1; l0<=1; l0++){
	  for(int l1=-1; l1<=1; l1++){
	    
	    std::vector<scalar_type> vec = vec_0;
	    
	    for(int d=0; d<DIMENSION; d++)
	      vec[d] += (l0*basis[0][d] + l1*basis[1][d]); 
	    
	    scalar_type distance = VECTOR_OPERATIONS::L2_NORM(vec, vec_1);
	    
	    if(distance < MIN_DISTANCE)
	      MIN_DISTANCE = distance;
	  }
	}
      }
      break;
      
    case 3:
      {
	for(int l0=-1; l0<=1; l0++){
	  for(int l1=-1; l1<=1; l1++){
	    for(int l2=-1; l2<=1; l2++){
	      
	      std::vector<scalar_type> vec = vec_0;
	      
	      for(int d=0; d<DIMENSION; d++)
		vec[d] += (l0*basis[0][d] + l1*basis[1][d] + l2*basis[2][d]); 
	      
	      scalar_type distance = VECTOR_OPERATIONS::L2_NORM(vec, vec_1);
	      
	      if(distance < MIN_DISTANCE)
		MIN_DISTANCE = distance;
	    }
	  }
	}
      }
      break;
      
    default:
      throw std::logic_error(__FUNCTION__);
    }

  return MIN_DISTANCE;
}

template<typename scalar_type>
bool cluster_operations::is_minimal(std::vector<scalar_type> R_vec,
				    std::vector<std::vector<scalar_type> >& basis)
{
  int DIMENSION = R_vec.size();
  
  std::vector<scalar_type> origin(DIMENSION,0.);

  scalar_type MIN_DISTANCE = minimal_distance(origin, R_vec, basis);
  
  bool minimal = VECTOR_OPERATIONS::L2_NORM(R_vec)>(MIN_DISTANCE+1.e-6)? false : true;
 
  return minimal;
}

template<typename scalar_type>
std::vector<std::vector<scalar_type> > cluster_operations::equivalent_vectors(std::vector<scalar_type> R_vec,
									      std::vector<std::vector<scalar_type> >& basis)
{
  const static scalar_type EPS          = 1.e-3;
  const static scalar_type ONE_PLUS_EPS = 1.+EPS;
  const static scalar_type ONE_MIN_EPS  = 1.-EPS;

  int DIMENSION = R_vec.size();

  R_vec = translate_inside_cluster(R_vec, basis);

  std::vector<scalar_type> origin(DIMENSION,0.);

  scalar_type MIN_DISTANCE = minimal_distance(origin, R_vec, basis);
  
  //bool IS_MINIMAL = L2_norm(R_vec)>(MIN_DISTANCE+1.e-6)? false : true;

  std::vector<std::vector<scalar_type> > r_min;
 
  switch(DIMENSION)
    {
    case 1:
      {
	for(int l0=-2; l0<=2; l0++){
	    	    
	  std::vector<scalar_type> vec = R_vec;
	  
	  for(int d=0; d<DIMENSION; d++)
	    vec[d] += l0*basis[0][d]; 
	  
	  scalar_type distance = sqrt(VECTOR_OPERATIONS::L2_NORM(vec));
	  
	  if(distance > sqrt(MIN_DISTANCE)*ONE_MIN_EPS -EPS &&
	     distance < sqrt(MIN_DISTANCE)*ONE_PLUS_EPS+EPS)
	    r_min.push_back(vec);
	}
      }	      
      break;

    case 2:
      {
	for(int l0=-2; l0<=2; l0++){
	  for(int l1=-2; l1<=2; l1++){
	    
	    std::vector<scalar_type> vec = R_vec;
	    
	    for(int d=0; d<DIMENSION; d++)
	      vec[d] += (l0*basis[0][d] + l1*basis[1][d]); 

	    scalar_type distance = sqrt(VECTOR_OPERATIONS::L2_NORM(vec));
	    
	    if(distance > sqrt(MIN_DISTANCE)*ONE_MIN_EPS -EPS && 
	       distance < sqrt(MIN_DISTANCE)*ONE_PLUS_EPS+EPS)
	      r_min.push_back(vec);
	  }
	}
      }     
      break;
      
    case 3:
      {
	for(int l0=-2; l0<=2; l0++){
	  for(int l1=-2; l1<=2; l1++){
	    for(int l2=-2; l2<=2; l2++){
	      
	      std::vector<scalar_type> vec = R_vec;
	      
	      for(int d=0; d<DIMENSION; d++)
		vec[d] += (l0*basis[0][d] + l1*basis[1][d] + l2*basis[2][d]); 
	      
	      scalar_type distance = sqrt(VECTOR_OPERATIONS::L2_NORM(vec));
	      
	      if(distance > sqrt(MIN_DISTANCE)*ONE_MIN_EPS -EPS && 
		 distance < sqrt(MIN_DISTANCE)*ONE_PLUS_EPS+EPS)
		r_min.push_back(vec);
	    }
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

  return r_min;
}

template<typename cluster_type, typename scalar_type>
std::pair<std::vector<scalar_type>, int> cluster_operations::find_closest_cluster_vector(const std::vector<scalar_type>& input_vec, const double tol)
{
  std::vector<std::vector<scalar_type>>& elements = cluster_type::get_elements();
  std::vector<std::vector<scalar_type>>& super_basis = cluster_type::get_super_basis_vectors();

  std::vector<scalar_type> input_vec_translated = translate_inside_cluster(input_vec, super_basis);

  if (elements.size() == 0) {
    throw std::logic_error(__FUNCTION__);
  }

  double min_distance = VECTOR_OPERATIONS::L2_NORM(input_vec_translated, elements[0]);
  int min_index = 0;
  
  for (int l=0; l < elements.size(); l++) {
    double distance = VECTOR_OPERATIONS::L2_NORM(input_vec_translated, elements[l]);
    if (distance < min_distance) {
      min_distance = distance;
      min_index = l;
    }
  }

  if (min_distance > tol) {
    throw std::logic_error(__FUNCTION__);
  }

  return std::make_pair(elements[min_index], min_index);
}


#endif
