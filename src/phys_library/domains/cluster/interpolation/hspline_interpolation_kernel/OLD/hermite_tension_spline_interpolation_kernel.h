//-*-C++-*-

#ifndef HERMITE_TENSION_SPLINE_INTERPOLATION_KERNEL_K_DMN_H
#define HERMITE_TENSION_SPLINE_INTERPOLATION_KERNEL_K_DMN_H

template<typename scalartype, typename source_k_dmn_t, typename target_k_dmn_t>
class hermite_tension_spline_interpolation_kernel
{};

/*!
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   This class implements a Hermite tension spline interpolation technique in momentum space.
 */
template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
class hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t> 
{
  typedef k_cluster<cluster_representation, base_cluster_type> k_cluster_type;

  const static int DIMENSION = k_cluster<cluster_representation, base_cluster_type>::DIMENSION;

  typedef k_cluster<cluster_representation, base_cluster_type> source_k_dmn_t;
  typedef r_cluster<cluster_representation, base_cluster_type> source_r_dmn_t;

  typedef dmn_0<source_k_dmn_t> k_dmn_t;

public:

  hermite_tension_spline_interpolation_kernel();
  ~hermite_tension_spline_interpolation_kernel();

  void reset();

  void execute(scalartype* input, scalartype* output);

private:

  void find_neighbours();

  void find_alphas(double* input);
  void find_alphas(std::complex<double>* input);

  double               evaluate_hermite_kernel(double x, double a);
  std::complex<double> evaluate_hermite_kernel(double x, std::complex<double> a);

  scalartype evaluate_hermite_kernel_at(int K_ind, std::vector<double>& k_vec);

  double               product(             double  x,              double  y);
  std::complex<double> product(std::complex<double> x, std::complex<double> y);

  void find_basis();
  void find_k_vecs();

private:

  typedef dmn<DIMENSION, int>     DIMENSION_domain;
  typedef dmn_0<DIMENSION_domain> DIMENSION_dmn;

  struct neighbours_domain
  {
  public:

    typedef neighbours_domain   this_type;
    typedef std::vector<double> element_type;

    static int& get_size(){
      static int size=0;
      return size;
    }

    static std::vector<element_type>& get_elements()
    {
      static std::vector<element_type> elements;
      return elements;
    }
  };

  typedef dmn_0<neighbours_domain> neighbours_dmn_t;

private:

  FUNC_LIB::function<int   , dmn_2<k_dmn_t, neighbours_dmn_t> > neighbours_index;

  FUNC_LIB::function<scalartype, dmn_2<k_dmn_t, neighbours_dmn_t> > strain;
  FUNC_LIB::function<scalartype, dmn_2<k_dmn_t, neighbours_dmn_t> > alphas;

  double k_basis    [DIMENSION*DIMENSION];
  double k_basis_inv[DIMENSION*DIMENSION];

  double*     k_vecs;
};

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::hermite_tension_spline_interpolation_kernel():
  k_vecs(NULL)
{
  find_neighbours();

  find_basis();

  find_k_vecs();
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::~hermite_tension_spline_interpolation_kernel()
{
  delete [] k_vecs;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::reset()
{
  delete [] k_vecs;

  find_neighbours();

  find_basis();

  find_k_vecs();
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::execute(scalartype* input, 
																	    scalartype* output)
{
  if(DIMENSION != 2)
    throw std::logic_error(__FUNCTION__);

  find_alphas(input);

  std::vector<double> K_vec (DIMENSION, 0.);

  std::vector<double> k_vec (DIMENSION, 0.);
  std::vector<double> k_diff(DIMENSION, 0.);

  for(int k_ind=0; k_ind<target_k_dmn_t::get_size(); ++k_ind){

    output[k_ind] = 0.;

    for(int d=0; d<DIMENSION; ++d)
      k_vec[d] = k_vecs[d+k_ind*DIMENSION];

//     VECTOR_OPERATIONS::PRINT(k_vec);
//     cout << endl << endl;
    
    for(int K_ind=0; K_ind<source_k_dmn_t::get_size(); ++K_ind){
      
      scalartype result = 0.;

      for(int l0=-1; l0<=1; ++l0){
	for(int l1=-1; l1<=1; ++l1){
       
	  K_vec = source_k_dmn_t::get_elements()[K_ind];
	  
	  for(int d=0; d<DIMENSION; ++d)
	    K_vec[d] += l0*k_cluster_type::get_basis()[0][d] + l1*k_cluster_type::get_basis()[1][d];
	  
	  k_diff = VECTOR_OPERATIONS::SUBTRACT(k_vec, K_vec);
	  
// 	  scalartype tmp = evaluate_hermite_kernel_at(K_ind, k_diff);
// 	  cout << "\t" << real(tmp) << "\t" << imag(tmp) << endl << endl;;

	  result += evaluate_hermite_kernel_at(K_ind, k_diff);
	}
      }

      // cout << "\t" << real(result) << "\t" << imag(result) << endl;

      output[k_ind] += product(result, input[K_ind]);
    }
  }
}




template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::find_alphas(double* input)
{
  const double A=4;

  alphas.reset();

  for(int K_ind=0; K_ind<source_k_dmn_t::get_size(); ++K_ind){
    for(int n_ind=0; n_ind<neighbours_dmn_t::dmn_size(); ++n_ind){
    
      int K_ind_0 = K_ind;
      int K_ind_1 = neighbours_index(K_ind_0, n_ind);
      int K_ind_2 = neighbours_index(K_ind_1, n_ind);

      double val = (input[K_ind_0]+input[K_ind_2])/2.-input[K_ind_1];
      val /= std::sqrt(VECTOR_OPERATIONS::L2_NORM(k_dmn_t::get_elements()[K_ind_0], 
						  k_dmn_t::get_elements()[K_ind_1]));

      alphas(K_ind, n_ind) = -1./(1.+std::exp(A*abs(val)));

//       cout << alphas(K_ind, n_ind) << "\t";
    }
//     cout << "\n";
  }
//   cout << "\n";
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::find_alphas(std::complex<double>* input)
{
  const double A=4;

  for(int K_ind=0; K_ind<source_k_dmn_t::get_size(); ++K_ind){
    for(int n_ind=0; n_ind<neighbours_dmn_t::dmn_size(); ++n_ind){
    
      int K_ind_0 = K_ind;
      int K_ind_1 = neighbours_index(K_ind_0, n_ind);
      int K_ind_2 = neighbours_index(K_ind_1, n_ind);

      std::complex<double> val = (input[K_ind_0]+input[K_ind_2])/2.-input[K_ind_1];
      
      real(alphas(K_ind, n_ind)) = -1./(1.+std::exp(A*abs(real(val))));
      imag(alphas(K_ind, n_ind)) = -1./(1.+std::exp(A*abs(imag(val))));
    }
  }
}



template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
double hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::evaluate_hermite_kernel(double x, 
																			      double a)
{
  double absX = fabs(x);

  if(absX <= 1.)
    return (a+2.)*pow(absX,3)-(a+3.)*pow(absX,2)+1.;

  if(absX > 1. and absX <= 2.)
    return a*pow(absX,3)-5.*a*pow(absX,2)+8.*a*absX-4.*a;

  return 0.;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
std::complex<double> hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::evaluate_hermite_kernel(double x, 
																					    std::complex<double> a)
{
  std::complex<double> result=0;

  double absX = fabs(x);

  if(absX <= 1.){
    real(result) = (real(a)+2.)*pow(absX,3)-(real(a)+3.)*pow(absX,2)+1.;
    imag(result) = (imag(a)+2.)*pow(absX,3)-(imag(a)+3.)*pow(absX,2)+1.;
  }

  if(absX > 1. and absX <= 2.){
    real(result) = real(a)*pow(absX,3)-5.*real(a)*pow(absX,2)+8.*real(a)*absX-4.*real(a);
    imag(result) = imag(a)*pow(absX,3)-5.*imag(a)*pow(absX,2)+8.*imag(a)*absX-4.*imag(a);
  }

  return result;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
scalartype hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::evaluate_hermite_kernel_at(int K_ind,
																				     std::vector<double>& k_vec)
{
  if(DIMENSION != 2)
    throw std::logic_error(__FUNCTION__);

  scalartype Hx, Hy;
  std::vector<double> x(DIMENSION,0);

  double     count  = 0.;
  scalartype result = 0.;

  for(int n0=0; n0<neighbours_dmn_t::dmn_size(); ++n0){
    for(int n1=0; n1<neighbours_dmn_t::dmn_size(); ++n1){

      std::vector<double> n0_vec = neighbours_dmn_t::get_elements()[n0];
      std::vector<double> n1_vec = neighbours_dmn_t::get_elements()[n1];

      if( (n0!=n1) and VECTOR_OPERATIONS::VOLUME(n0_vec, n1_vec)>1.e-6){

	VECTOR_OPERATIONS::COORDINATES(n0_vec, n1_vec, k_vec, x);
	
	if(x[0]>=0. && x[1]>=0. && x[0]<=2. && x[1]<=2.){

	  Hx = evaluate_hermite_kernel(x[0], alphas(K_ind, n0));
	  Hy = evaluate_hermite_kernel(x[1], alphas(K_ind, n1));

// 	  cout << alphas(K_ind, n0) << "\t" << Hx << "\t" << Hy << endl;
	  count  += 1.;
	  result += product(Hx,Hy);
	}
      }
    }
  }

//   if(count>0)
//     cout << "\t" << real(result/count) << "\t" << imag(result/count) << endl;

  if(count>0)
    return result/count;
  else
    return 0.;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
inline double hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::product(double x, double y)
{
  return x*y;
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
inline std::complex<double> hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::product(std::complex<double> x,
																				   std::complex<double> y)
{
  return std::complex<double>(real(x)*real(y), imag(x)*imag(y));
}


template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::find_neighbours()
{
  tetrahedron_mesh<k_cluster_type> tet_mesh(1);
  
  neighbours_domain::get_size() = 0;
  neighbours_domain::get_elements().resize(0);

  for(size_t f_ind=0; f_ind<tet_mesh.get_facets().size(); ++f_ind){
    
    std::vector<double> k(DIMENSION,0.);
    
    for(size_t k_ind=0; k_ind<tet_mesh.get_facets()[f_ind].index.size(); ++k_ind)
      k = VECTOR_OPERATIONS::ADD(k, tet_mesh.get_simplices()[tet_mesh.get_facets()[f_ind].index[k_ind]].k_vec);

    neighbours_domain::get_size() += 1;
    neighbours_domain::get_elements().push_back(k);
  }

  assert(neighbours_domain::get_elements().size() == tet_mesh.get_facets().size());

  neighbours_index.reset();
  alphas.reset();

  for(int K_ind=0; K_ind<source_k_dmn_t::get_size(); ++K_ind){
    for(int n_ind=0; n_ind<neighbours_dmn_t::dmn_size(); ++n_ind){
      
      std::vector<double> k(DIMENSION,0.);

      k = VECTOR_OPERATIONS::ADD(source_k_dmn_t  ::get_elements()[K_ind], 
				 neighbours_dmn_t::get_elements()[n_ind]);

      k = source_k_dmn_t::back_inside_cluster(k);

      int index=-1;
      for(int l=0; l<source_k_dmn_t::get_size(); ++l){
	if(VECTOR_OPERATIONS::L2_NORM(k, source_k_dmn_t::get_elements()[l])){
	  index=l;
	  break;
	}
      }

      assert(index>-1);
      neighbours_index(K_ind, n_ind) = index;
    }
  }
}

template<typename scalartype, cluster_representation_type cluster_representation, typename base_cluster_type, typename target_k_dmn_t>
void hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::find_basis()
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
void hermite_tension_spline_interpolation_kernel<scalartype, k_cluster<cluster_representation, base_cluster_type>, target_k_dmn_t>::find_k_vecs()
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


#endif

























