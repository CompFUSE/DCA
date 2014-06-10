//-*-C++-*-

#ifndef HERMITE_SPLINE_INTERPOLATION_KERNEL_K_DMN_2ND_ORDER_H
#define HERMITE_SPLINE_INTERPOLATION_KERNEL_K_DMN_2ND_ORDER_H

template<typename scalartype, typename source_k_dmn_t, typename target_k_dmn_t>
class hspline_interpolation_kernel_second_order
{};

/*!
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   This class implements a Hermite spline interpolation technique in momentum space.
 */
template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
class hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t> 
{
  typedef k_cluster<cluster_representation, base_cluster_type> k_cluster_type;

  const static int DIMENSION = k_cluster<cluster_representation, base_cluster_type>::DIMENSION;

  typedef k_cluster<cluster_representation, base_cluster_type> source_k_dmn_t;
  typedef r_cluster<cluster_representation, base_cluster_type> source_r_dmn_t;

public:

  hspline_interpolation_kernel_second_order();
  hspline_interpolation_kernel_second_order(int A);

  ~hspline_interpolation_kernel_second_order();

  void reset();

  scalartype* get_interpolation_matrix();

  void execute(scalartype* input, scalartype* output);

  void execute(scalartype* input, scalartype* output, int n);

private:

//   void find_neighbours();

  void find_basis();
  void find_k_vecs();

  double evaluate_hermite_kernel(double x);
  double evaluate_hermite_kernel_at(std::vector<double>& k_vec);

  void construct_interpolation_matrix();

private:

  double a;

//   std::vector<std::vector<double> > neighbours;

  double k_basis    [DIMENSION*DIMENSION];
  double k_basis_inv[DIMENSION*DIMENSION];

  double*     k_vecs;

//   scalartype* K_K_matrix;
  scalartype* interpolation_matrix;
};

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::hspline_interpolation_kernel_second_order():
  a(-0.5),

//   neighbours(0),

  k_vecs(NULL),
  interpolation_matrix(NULL)
{
//   find_neighbours();

  find_basis();

  find_k_vecs();

  construct_interpolation_matrix();
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::hspline_interpolation_kernel_second_order(int A):
  a(A),

//   neighbours(0),

  k_vecs(NULL),
  interpolation_matrix(NULL)
{
//   find_neighbours();

  find_basis();

  find_k_vecs();

  construct_interpolation_matrix();
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::~hspline_interpolation_kernel_second_order()
{
  delete [] k_vecs;
  delete [] interpolation_matrix;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::reset()
{
  delete [] k_vecs;
  delete [] interpolation_matrix;

//   find_neighbours();

  find_basis();

  find_k_vecs();

  construct_interpolation_matrix();
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
scalartype* hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::get_interpolation_matrix()
{
  return interpolation_matrix;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::execute(scalartype* input, 
																	  scalartype* output)
{
  gemm_plan<scalartype> gemm(target_k_dmn_t::get_size(), source_k_dmn_t::get_size(), 1.);
    
  gemm.A = interpolation_matrix;
  gemm.B = input;
  gemm.C = output;
    
  gemm.execute_plan();
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::execute(scalartype* input, 
																	  scalartype* output,
																	  int n)
{
  gemm_plan<scalartype> gemm(target_k_dmn_t::get_size(), source_k_dmn_t::get_size(), n);
  
  gemm.A = interpolation_matrix;
  gemm.B = input;
  gemm.C = output;

  gemm.execute_plan();
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::find_basis()
{
  for(int d0=0; d0<DIMENSION; ++d0)
    for(int d1=0; d1<DIMENSION; ++d1)
      k_basis[d0+d1*DIMENSION] = source_k_dmn_t::get_basis()[d0][d1];

  {
    invert_plan<double> invert_pln(DIMENSION);
    memcpy(invert_pln.Matrix, k_basis, sizeof(double)*DIMENSION*DIMENSION);
    invert_pln.execute_plan();
    memcpy(k_basis_inv, invert_pln.inverted_matrix, sizeof(double)*DIMENSION*DIMENSION);
  }

  for(int d0=0; d0<DIMENSION; ++d0)
    for(int d1=0; d1<DIMENSION; ++d1)
      k_basis[d1+d0*DIMENSION] = source_k_dmn_t::get_basis()[d0][d1];
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::find_k_vecs()
{
  double* k_vecs_target = new double[DIMENSION*target_k_dmn_t::get_size()];
  double* k_vecs_source = new double[DIMENSION*target_k_dmn_t::get_size()];

  for(int l=0; l<target_k_dmn_t::get_size(); ++l)
    for(int d=0; d<DIMENSION; ++d)
      k_vecs_target[d+l*DIMENSION] = target_k_dmn_t::get_elements()[l][d];

  {
    gemm_plan<double> gemm(DIMENSION, DIMENSION, target_k_dmn_t::get_size());
    
    gemm.A = k_basis_inv;
    gemm.B = k_vecs_target;
    gemm.C = k_vecs_source;

    gemm.execute_plan();
  }

  for(int l=0; l<target_k_dmn_t::get_size(); ++l){
    for(int d=0; d<DIMENSION; ++d){

      while(k_vecs_source[d+l*DIMENSION] > 1.-1.e-6)
	k_vecs_source[d+l*DIMENSION] -= 1.;

      while(k_vecs_source[d+l*DIMENSION] < -1.e-6)
	k_vecs_source[d+l*DIMENSION] += 1.;

      assert(k_vecs_source[d+l*DIMENSION] > 0.-1.e-6);
      assert(k_vecs_source[d+l*DIMENSION] < 1.+1.e-6);
    }
  }

  k_vecs = new double[DIMENSION*target_k_dmn_t::get_size()];

  {
    gemm_plan<double> gemm(DIMENSION, DIMENSION, target_k_dmn_t::get_size());
    
    gemm.A = k_basis;
    gemm.B = k_vecs_source;
    gemm.C = k_vecs;

    gemm.execute_plan();
  }

  delete [] k_vecs_source;
  delete [] k_vecs_target;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
double hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::evaluate_hermite_kernel(double x)
{
  double absX = fabs(x);

  if(absX <= 1.)
    return (a+2)*pow(absX,3)-(a+3)*pow(absX,2)+1;

  if(absX > 1. and absX <= 2.)
    return a*pow(absX,3)-5.*a*pow(absX,2)+8*a*absX-4*a;

  return 0.;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
double hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::evaluate_hermite_kernel_at(std::vector<double>& k_vec)
{
  typedef tetrahedron_neighbour_domain<k_cluster_type> tet_dmn_t;

  std::vector<std::vector<double> >& k_vecs = tet_dmn_t::get_elements();

  std::vector<double> x(DIMENSION,0);

//   double count  = 0.;
  double result = 0.;

  for(size_t f_ind=0; f_ind<tet_dmn_t::get_facets().size(); ++f_ind){

    switch(DIMENSION)
      {
      case 2:
	{
	  std::vector<double> k0 = k_vecs[tet_dmn_t::get_facets()[f_ind][0]];
	  std::vector<double> k1 = k_vecs[tet_dmn_t::get_facets()[f_ind][1]];

	  assert(VECTOR_OPERATIONS::VOLUME(k0, k1)>1.e-6);
	  
	  VECTOR_OPERATIONS::COORDINATES(k0, k1, k_vec, x);
	  
	  bool lies_in_positive_part=true;
	  for(int d=0; d<DIMENSION; d++)
	    if(x[d]<-1.e-6)
	      lies_in_positive_part=false;
	  
	  if(lies_in_positive_part)
	    {
	      std::vector<double> k(DIMENSION,0);
	      for(int d=0; d<DIMENSION; d++)
		k[d] = k0[d]+k1[d];

	      double phi = VECTOR_OPERATIONS::DOT_PRODUCT(k, k_vec)/VECTOR_OPERATIONS::DOT_PRODUCT(k, k); 

	      assert(phi>-1.e-6);
	      result = evaluate_hermite_kernel(phi);
	    }
	}
	break;

      default:
	throw std::logic_error(__FUNCTION__);
      }
  }

  return result;
//   if(count>0)
//     return result/count;
//   else
//     return 0.;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::construct_interpolation_matrix()
{
  if(DIMENSION != 2)
    throw std::logic_error(__FUNCTION__);

  interpolation_matrix = new scalartype[target_k_dmn_t::get_size()*source_k_dmn_t::get_size()];

  std::vector<double> K_vec (DIMENSION, 0.);

  std::vector<double> k_vec (DIMENSION, 0.);
  std::vector<double> k_diff(DIMENSION, 0.);

  for(int k_ind=0; k_ind<target_k_dmn_t::get_size(); ++k_ind){

    for(int d=0; d<DIMENSION; ++d)
      k_vec[d] = k_vecs[d+k_ind*DIMENSION];

    for(int K_ind=0; K_ind<source_k_dmn_t::get_size(); ++K_ind){

      scalartype result = 0.;

      for(int l0=-1; l0<=1; ++l0){
	for(int l1=-1; l1<=1; ++l1){
       
	  K_vec = source_k_dmn_t::get_elements()[K_ind];
	  
	  for(int d=0; d<DIMENSION; ++d)
	    K_vec[d] += l0*k_cluster_type::get_basis()[0][d] + l1*k_cluster_type::get_basis()[1][d];
	  
	  k_diff = VECTOR_OPERATIONS::SUBTRACT(k_vec, K_vec);
	  
	  result += evaluate_hermite_kernel_at(k_diff);
	}
      }

      interpolation_matrix[k_ind+target_k_dmn_t::get_size()*K_ind] = result;
    }
  }

  for(int k_ind=0; k_ind<target_k_dmn_t::get_size(); ++k_ind){
    for(int K_ind=0; K_ind<source_k_dmn_t::get_size(); ++K_ind){
      cout << real(interpolation_matrix[k_ind+target_k_dmn_t::get_size()*K_ind]) << "\t";
    }
    cout << endl;
  }
  cout << endl;
}


#endif


/*
template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
double hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::evaluate_hermite_kernel_at(std::vector<double>& k_vec)
{
  if(DIMENSION != 2)
    throw std::logic_error(__FUNCTION__);

  double Hx, Hy;
  std::vector<double> x(DIMENSION,0);

  double count  = 0.;
  double result = 0.;

  for(size_t n0=0; n0<neighbours.size(); ++n0){
    for(size_t n1=0; n1<neighbours.size(); ++n1){

      if(VECTOR_OPERATIONS::VOLUME(neighbours[n0], neighbours[n1])>1.e-6){

	VECTOR_OPERATIONS::COORDINATES(neighbours[n0], neighbours[n1], k_vec, x);
	
 	Hx = evaluate_hermite_kernel(x[0]);
 	Hy = evaluate_hermite_kernel(x[1]);

	count  += 1.;
	result += Hx*Hy;
      }
    }
  }

  return result/count;
}
 */

/*
template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
double hspline_interpolation_kernel_second_order<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::evaluate_hermite_kernel_at(std::vector<double>& k_vec)
{
  double dot_prod, norm_neigh;

  double result=1.;

  double x=0.;

  for(size_t k_ind=0; k_ind<neighbours.size(); ++k_ind){
    
    norm_neigh = VECTOR_OPERATIONS::L2_NORM(neighbours[k_ind]);
    
    dot_prod = VECTOR_OPERATIONS::DOT_PRODUCT(k_vec, neighbours[k_ind]);
	
    x = dot_prod/norm_neigh;
    
    if(x >= 0.){
	  
      if(x <= 1.)
	result *= (a+2)*pow(x,3)-(a+3)*pow(x,2)+1;
      else
	{
	  if(x > 1. and x <= 2)
	    result *= a*pow(x,3)-5.*a*pow(x,2)+8*a*x-4*a;
	  else
	    result *= 0.;
	}
    }
  }

  return result;
}
*/
