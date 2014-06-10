//-*-C++-*-

#ifndef HERMITE_SPLINE_INTERPOLATION_KERNEL_K_DMN_H
#define HERMITE_SPLINE_INTERPOLATION_KERNEL_K_DMN_H

/*!
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   This class implements a Hermite spline interpolation technique in momentum space.
 */
template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
class hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t> 
{
  typedef k_cluster<cluster_representation, base_cluster_type> k_cluster_type;

  const static int DIMENSION = k_cluster<cluster_representation, base_cluster_type>::DIMENSION;

  typedef k_cluster<cluster_representation, base_cluster_type> source_k_dmn_t;
  typedef r_cluster<cluster_representation, base_cluster_type> source_r_dmn_t;

public:

  hspline_interpolation_kernel();
  hspline_interpolation_kernel(double A);

  ~hspline_interpolation_kernel();

  void reset();

  scalartype* get_interpolation_matrix();

  void execute(scalartype* input, scalartype* output);

  void execute(scalartype* input, scalartype* output, int n);

  void execute_on_transpose(scalartype* input, scalartype* output, int n);

private:

  void find_neighbours();

  void find_basis();
  void find_k_vecs();

  double evaluate_hermite_kernel(double x);
  double evaluate_hermite_kernel_at(std::vector<double>& k_vec);

  void construct_interpolation_matrix();

  static double volume(double* b0, double* b1);
  static double volume(double* b0, double* b1, double* b);

  static void coordinates(double* b0, double* b1,             double* r, double* x); 
  static void coordinates(double* b0, double* b1, double* b2, double* r, double* x); 

private:

  double a;

  std::vector<std::vector<double> > neighbours;

  double k_basis    [DIMENSION*DIMENSION];
  double k_basis_inv[DIMENSION*DIMENSION];

  double*     k_vecs;

  scalartype* interpolation_matrix;
};

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::hspline_interpolation_kernel():
  a(-0.5),

  neighbours(0),

  k_vecs(NULL),
  interpolation_matrix(NULL)
{
  find_neighbours();

  find_basis();

  find_k_vecs();

  construct_interpolation_matrix();
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::hspline_interpolation_kernel(double A):
  a(A),

  neighbours(0),

  k_vecs(NULL),
  interpolation_matrix(NULL)
{
  find_neighbours();

  find_basis();

  find_k_vecs();

  construct_interpolation_matrix();
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::~hspline_interpolation_kernel()
{
  delete [] k_vecs;
  delete [] interpolation_matrix;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::reset()
{
  delete [] k_vecs;
  delete [] interpolation_matrix;

  find_k_vecs();

  construct_interpolation_matrix();
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
scalartype* hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::get_interpolation_matrix()
{
  return interpolation_matrix;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::execute(scalartype* input, 
															     scalartype* output)
{
  gemm_plan<scalartype> gemm(target_k_dmn_t::get_size(), source_k_dmn_t::get_size(), 1.);
    
  gemm.A = interpolation_matrix;
  gemm.B = input;
  gemm.C = output;
    
  gemm.execute_plan();
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::execute(scalartype* input, 
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
void hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::execute_on_transpose(scalartype* input, scalartype* output, int n)
{
  gemm_plan<scalartype> gemm(n, source_k_dmn_t::get_size(), target_k_dmn_t::get_size());
  
  gemm.TRANSA = 'N';
  gemm.TRANSB = 'T';

  gemm.LDA = n;
  gemm.LDB = target_k_dmn_t::get_size();
  gemm.LDC = n;

  gemm.A = input;
  gemm.B = interpolation_matrix;
  gemm.C = output;

  gemm.execute_plan();
}





template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::find_neighbours()
{
  tetrahedron_mesh<k_cluster_type> tet_mesh(1);
  
  neighbours.resize(0);

  for(size_t f_ind=0; f_ind<tet_mesh.get_facets().size(); ++f_ind){
    
    std::vector<double> k(DIMENSION,0.);
    
    for(size_t k_ind=0; k_ind<tet_mesh.get_facets()[f_ind].index.size(); ++k_ind)
      k = VECTOR_OPERATIONS::ADD(k, tet_mesh.get_simplices()[tet_mesh.get_facets()[f_ind].index[k_ind]].k_vec);
    
    neighbours.push_back(k);
  }

  assert(neighbours.size() == tet_mesh.get_facets().size());
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::find_basis()
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
void hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::find_k_vecs()
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
inline double hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::volume(double* b0, double* b1)
{
  return (b0[0]*b1[1]-b0[1]*b1[0]);
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
inline double hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::volume(double* b0, double* b1, double* b2)
{
  return (-b0[2]*b1[1]*b2[0] + b0[1]*b1[2]*b2[0] + b0[2]*b1[0]*b2[1] - b0[0]*b1[2]*b2[1] - b0[1]*b1[0]*b2[2] + b0[0]*b1[1]*b2[2]);
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
inline void hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::coordinates(double* b0, double* b1, 
																	double* r, double* x)
{
  double inv_det = 1./(b0[0]*b1[1]-b0[1]*b1[0]);

  x[0] =  b1[1]*r[0] - b1[0]*r[1];
  x[1] = -b0[1]*r[0] + b0[0]*r[1];

  x[0] *= inv_det;
  x[1] *= inv_det;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
inline void hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::coordinates(double* b0, double* b1, double* b2, 
																	double* r, double* x)
{
  double inv_det = 1./(-b0[2]*b1[1]*b2[0] + b0[1]*b1[2]*b2[0] + b0[2]*b1[0]*b2[1] - b0[0]*b1[2]*b2[1] - b0[1]*b1[0]*b2[2] + b0[0]*b1[1]*b2[2]);

  x[0] = (-b1[2]*b2[1] + b1[1]*b2[2])*r[0] + ( b1[2]*b2[0] - b1[0]*b2[2])*r[1] + (-b1[1]*b2[0] + b1[0]*b2[1])*r[2];
  x[1] = ( b0[2]*b2[1] - b0[1]*b2[2])*r[0] + (-b0[2]*b2[0] + b0[0]*b2[2])*r[1] + ( b0[1]*b2[0] - b0[0]*b2[1])*r[2];
  x[2] = (-b0[2]*b1[1] + b0[1]*b1[2])*r[0] + ( b0[2]*b1[0] - b0[0]*b1[2])*r[1] + (-b0[1]*b1[0] + b0[0]*b1[1])*r[2];

  x[0] *= inv_det;
  x[1] *= inv_det;
  x[2] *= inv_det;
}


template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
inline double hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::evaluate_hermite_kernel(double x)
{
  double absX = fabs(x);

  if(absX>2)
    return 0.;

  if(absX>1)
    return a*pow(absX,3)-5.*a*pow(absX,2)+8.*a*absX-4.*a;

  return (a+2.)*pow(absX,3)-(a+3.)*pow(absX,2)+1.;

  /*
  if(absX <= 1.)
    return (a+2)*pow(absX,3)-(a+3)*pow(absX,2)+1;

  if(absX > 1. and absX <= 2.)
    return a*pow(absX,3)-5.*a*pow(absX,2)+8*a*absX-4*a;

  return 0.;
  */
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
inline double hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::evaluate_hermite_kernel_at(std::vector<double>& k_vec)
{
  std::vector<double> x(DIMENSION,0);

  double count  = 0.;
  double result = 0.;

  switch(DIMENSION)
    {
    case 1:
      {
	double Hx;
	for(size_t n0=0; n0<neighbours.size(); ++n0){
	    
	  x[0] = k_vec[0]/fabs(neighbours[n0][0]);
  
	  Hx = evaluate_hermite_kernel(x[0]);
	  
	  count  += 1.;
	  result += Hx;
	}
      }
      break;

    case 2:
      {
	double Hx, Hy;
	for(size_t n0=0; n0<neighbours.size(); ++n0){
	  for(size_t n1=0; n1<neighbours.size(); ++n1){
	    
	    if(fabs(volume(&neighbours[n0][0], &neighbours[n1][0])) > 1.e-6){

	      coordinates(&neighbours[n0][0], &neighbours[n1][0], &k_vec[0], &x[0]);
	      
	      Hx = evaluate_hermite_kernel(x[0]);
	      Hy = evaluate_hermite_kernel(x[1]);
	      
	      count  += 1.;
	      result += Hx*Hy;	    
	    }
	  }
	}
      }
      break;

    case 3:
      {
	double Hx, Hy, Hz;
	for(size_t n0=0; n0<neighbours.size(); ++n0){
	  for(size_t n1=0; n1<neighbours.size(); ++n1){
	    for(size_t n2=0; n2<neighbours.size(); ++n2){
	    
	      if(fabs(volume(&neighbours[n0][0], &neighbours[n1][0], &neighbours[n2][0])) > 1.e-6){
// 	      if(VECTOR_OPERATIONS::VOLUME(neighbours[n0], neighbours[n1], neighbours[n2])>1.e-6){

		coordinates(&neighbours[n0][0], &neighbours[n1][0], &neighbours[n2][0], &k_vec[0], &x[0]);
// 		VECTOR_OPERATIONS::COORDINATES(neighbours[n0], neighbours[n1], neighbours[n2], k_vec, x);
		
		Hx = evaluate_hermite_kernel(x[0]);
		Hy = evaluate_hermite_kernel(x[1]);
		Hz = evaluate_hermite_kernel(x[2]);
		
		count  += 1.;
		result += Hx*Hy*Hz;
	      }
	    }
	  }
	}
      }
      break;
      
    default :
      throw std::logic_error(__FUNCTION__);
    };

  if(count>0)
    return result/count;
  else
    return 0.;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hspline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::construct_interpolation_matrix()
{
  int N = extended_k_domain<source_k_dmn_t>::get_size();

  interpolation_matrix = new scalartype[target_k_dmn_t::get_size()*source_k_dmn_t::get_size()];

  for(int i=0; i<target_k_dmn_t::get_size()*source_k_dmn_t::get_size(); ++i)
    interpolation_matrix[i] = 0.;

  std::vector<double> K_vec (DIMENSION, 0.);

  std::vector<double> k_vec (DIMENSION, 0.);
  std::vector<double> k_diff(DIMENSION, 0.);

  for(int k_ind=0; k_ind<target_k_dmn_t::get_size(); ++k_ind){

    for(int d=0; d<DIMENSION; ++d)
      k_vec[d] = k_vecs[d+k_ind*DIMENSION];
        
    for(int n_ind=0; n_ind<N; ++n_ind){
      
      K_vec = extended_k_domain<source_k_dmn_t>::get_elements()[n_ind];
      
      int K_ind = extended_k_domain<source_k_dmn_t>::get_cluster_index(n_ind);

      k_diff = VECTOR_OPERATIONS::SUBTRACT(k_vec, K_vec);

      interpolation_matrix[k_ind+target_k_dmn_t::get_size()*K_ind] += evaluate_hermite_kernel_at(k_diff);
    }
  }
}

#endif





















