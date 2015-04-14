//-*-C++-*-

#ifndef WANNIER_INTERPOLATION_KERNEL_K_TO_k_NFFT_NEW_H
#define WANNIER_INTERPOLATION_KERNEL_K_TO_k_NFFT_NEW_H

/*! \class   wannier_interpolation_kernel
 *  \ingroup INTERPOLATION
 *
 *  \author  Peter Staar
 *  \brief   This class implements a Wannier interpolation technique, using the nfft library.
 */
//template<typename source_base_cluster_type, typename target_dmn_type>
template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
class wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>
{
  //#include "type_definitions.h"

//   const static int DIMENSION = source_base_cluster_type::DIMENSION;

//   typedef r_cluster<FULL, source_base_cluster_type> source_r_cluster_type;
//   typedef k_cluster<FULL, source_base_cluster_type> source_k_cluster_type;

//   typedef dmn_0<source_r_cluster_type> source_r_dmn_t;
//   typedef dmn_0<source_k_cluster_type> source_k_dmn_t;

  typedef cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S> k_cluster_type;
  
  const static int DIMENSION = k_cluster_type::DIMENSION;

  typedef typename k_cluster_type::dual_type source_r_cluster_type;
  typedef typename k_cluster_type::this_type source_k_cluster_type;

  typedef dmn_0<source_r_cluster_type> source_r_dmn_t;
  typedef dmn_0<source_k_cluster_type> source_k_dmn_t;

  typedef dmn_0<target_dmn_type> target_k_dmn_t;

public:

  //typedef dmn_0<centered_r_LDA> centered_r_LDA_dmn_t;

public:

  wannier_interpolation_kernel();
  ~wannier_interpolation_kernel();

  void reset_input();
  void reset_output();

  void set(std::complex<double>* input_ptr);
  void get(std::complex<double>* output_ptr);

  void execute(std::complex<double>* input_ptr,
	       std::complex<double>* output_ptr);

  std::complex<double>& get_F_r(int i);

private:

  void check_grid_sizes();

  void initialize_centered_r_cluster();

  void initialize_nfft_K_2_R();
  void initialize_nfft_R_2_k();

  void initialize_cut_off();

  void FT_to_centered_function_NFFT();

  void FT_F_K__to__F_R(std::complex<double>* input_ptr);

  void FT_F_R__to__F_k(std::complex<double>* output_ptr);

private:

  struct centered_r_cluster {

    typedef centered_r_cluster      this_type;
    typedef std::vector<double> element_type;

    static int get_size() {
      return get_elements().size();
      //return r_LDA::dmn_size();
    }

    static std::vector<element_type>& get_elements() {
      static std::vector<element_type> elements(0);
      return elements;
    }
  };

public:

  typedef dmn_0<centered_r_cluster> centered_r_cluster_dmn_t;

private:

  static bool INITIALIZED;

  static FUNC_LIB::function<int, centered_r_cluster_dmn_t> lies_within_cutoff; 

  std::vector<int> grid_size;

  nfft_plan nfft_K_2_R;
  nfft_plan nfft_R_2_k;
  
  FUNC_LIB::function<std::complex<double>, centered_r_cluster_dmn_t> F_R;
};

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
bool wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::INITIALIZED = false;

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
FUNC_LIB::function<int, typename wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::centered_r_cluster_dmn_t> 
wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::lies_within_cutoff("cutoff");

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::wannier_interpolation_kernel():
  //grid_size(LDA_cluster_type::grid_size()),
  //grid_size(k_LDA::parameter_type::get_dimensions()),

  grid_size(DIMENSION, 32),

  nfft_K_2_R(),
  nfft_R_2_k(),

  F_R("wannier_interpolation_kernel__F_r")
{
  for(int i=0; i<DIMENSION; ++i)
    grid_size[i] = grid_size[i]<=4 ? 6 : grid_size[i];

  check_grid_sizes();

  if(!INITIALIZED)
    {
      initialize_centered_r_cluster();
      
      F_R               .reset();
      lies_within_cutoff.reset();
      
      initialize_cut_off();
      
      INITIALIZED = true;
    }

  initialize_nfft_K_2_R();
  initialize_nfft_R_2_k();
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::~wannier_interpolation_kernel()
{
  nfft_finalize(&nfft_K_2_R);
  nfft_finalize(&nfft_R_2_k);
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::check_grid_sizes()
{
  int size = 1;
  for(int i=0; i<DIMENSION; ++i)
    size *= grid_size[i];

  if(size < source_k_dmn_t::dmn_size()){
    cout << "\n\n\t INCREASE LDA-grid FOR WANNIER-INTERPOLATION!!! \n\n\n";
    throw std::logic_error(__FUNCTION__);
  }
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
std::complex<double>& wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::get_F_r(int i)
{
  return F_R(i);
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::reset_input()
{
  nfft_finalize(&nfft_K_2_R);
  initialize_nfft_K_2_R();
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::reset_output()
{
  nfft_finalize(&nfft_R_2_k);
  initialize_nfft_R_2_k();
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::execute(std::complex<double>* input_ptr,
												       std::complex<double>* output_ptr)
{
  FT_F_K__to__F_R(input_ptr);

  FT_F_R__to__F_k(output_ptr);
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::set(std::complex<double>* input_ptr)
{
  FT_F_K__to__F_R(input_ptr);
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::get(std::complex<double>* output_ptr)
{
  FT_F_R__to__F_k(output_ptr);
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::initialize_nfft_K_2_R()
{
  nfft_init(&nfft_K_2_R, DIMENSION, &grid_size[0], source_k_dmn_t::dmn_size());

  std::vector<std::vector<double> >& collection_k_vecs = source_k_dmn_t::get_elements();

  std::vector<std::vector<double> > collection_k_vecs_affine(source_k_dmn_t::dmn_size(), std::vector<double>(DIMENSION, 0));

  for(int j=0; j<source_k_dmn_t::dmn_size(); j++){
    //collection_k_vecs_affine[j] = source_k_cluster_type::get_affine_coordinate(collection_k_vecs[j]);
    collection_k_vecs_affine[j] = VECTOR_OPERATIONS::COORDINATES(collection_k_vecs[j], 
								 source_k_cluster_type::get_super_basis_vectors());

    //VECTOR_OPERATIONS::PRINT(collection_k_vecs_affine[j]); cout<<endl;
  }
  //cout<<endl;

  for(int j=0; j<source_k_dmn_t::dmn_size(); j++){
    for(int i=0; i<DIMENSION; i++){

      while(collection_k_vecs_affine[j][i] < -1./2.)
	collection_k_vecs_affine[j][i] += 1.;

      while(collection_k_vecs_affine[j][i] > 1./2.-1.e-6)
	collection_k_vecs_affine[j][i] -= 1.;
    }
  }
  
  for(int j=0; j<source_k_dmn_t::dmn_size(); j++){
    for(int i=0; i<DIMENSION; i++){
      nfft_K_2_R.x[j*DIMENSION + i] = collection_k_vecs_affine[j][i];
	
      assert(nfft_K_2_R.x[j*DIMENSION + i] >= -0.5-1.e-6);
      assert(nfft_K_2_R.x[j*DIMENSION + i] <   0.5);
    }
  }

  if(nfft_K_2_R.nfft_flags & PRE_ONE_PSI) 
    nfft_precompute_one_psi(&nfft_K_2_R);
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::initialize_nfft_R_2_k()
{
  nfft_init(&nfft_R_2_k, DIMENSION, &grid_size[0], target_k_dmn_t::dmn_size());

  std::vector<std::vector<double> >& collection_k_vecs = target_k_dmn_t::get_elements();

  std::vector<std::vector<double> > collection_k_vecs_affine(target_k_dmn_t::dmn_size(), std::vector<double>(DIMENSION, 0));

  for(int j=0; j<target_k_dmn_t::dmn_size(); j++){
    //collection_k_vecs_affine[j] = source_k_cluster_type::get_affine_coordinate(collection_k_vecs[j]);
    collection_k_vecs_affine[j] = VECTOR_OPERATIONS::COORDINATES(collection_k_vecs[j],
								 source_k_cluster_type::get_super_basis_vectors());

    //VECTOR_OPERATIONS::PRINT(collection_k_vecs_affine[j]); cout<<endl;
  }
  //cout<<endl;

  for(int j=0; j<target_k_dmn_t::dmn_size(); j++){
    for(int i=0; i<DIMENSION; i++){

      while(collection_k_vecs_affine[j][i] < -1./2.)
	collection_k_vecs_affine[j][i] += 1.;

      while(collection_k_vecs_affine[j][i] > 1./2.-1.e-6)
	collection_k_vecs_affine[j][i] -= 1.;
    }
  }
  
  for(int j=0; j<target_k_dmn_t::dmn_size(); j++){
    for(int i=0; i<DIMENSION; i++){
      nfft_R_2_k.x[j*DIMENSION + i] = collection_k_vecs_affine[j][i];
	
      assert(nfft_R_2_k.x[j*DIMENSION + i] >= -0.5-1.e-6);
      assert(nfft_R_2_k.x[j*DIMENSION + i] <   0.5);
    }
  }

  if(nfft_R_2_k.nfft_flags & PRE_ONE_PSI) 
    nfft_precompute_one_psi(&nfft_R_2_k);
}




template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::FT_F_K__to__F_R(std::complex<double>* input_ptr)
{
  assert(source_k_dmn_t::dmn_size()>0);

  for(int K_ind=0; K_ind<source_k_dmn_t::dmn_size(); K_ind++){

//     VECTOR_OPERATIONS::PRINT(source_k_dmn_t::get_elements()[K_ind]);
//     cout << real(input_ptr[K_ind]) << "\t" << imag(input_ptr[K_ind]) << endl;

    nfft_K_2_R.f[K_ind][0] = real(input_ptr[K_ind]);
    nfft_K_2_R.f[K_ind][1] = imag(input_ptr[K_ind]);
  }
//   cout << endl;

  nfft_adjoint(&nfft_K_2_R);
  
  for(int R_ind=0; R_ind<centered_r_cluster_dmn_t::dmn_size(); R_ind++){

    if(lies_within_cutoff(R_ind)>0){
	F_R(R_ind).real( nfft_K_2_R.f_hat[R_ind][0]/double(source_k_dmn_t::dmn_size()*lies_within_cutoff(R_ind)) );
	F_R(R_ind).imag( nfft_K_2_R.f_hat[R_ind][1]/double(source_k_dmn_t::dmn_size()*lies_within_cutoff(R_ind)) ); 
    }
    else 
      F_R(R_ind) = 0.;

//     VECTOR_OPERATIONS::PRINT(centered_r_cluster_dmn_t::get_elements()[R_ind]);
//     cout << real(F_R(R_ind)) << "\t" << imag(F_R(R_ind)) << endl;
  }
//   cout << endl;
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::FT_F_R__to__F_k(std::complex<double>* output_ptr)
{
  for(int r_ind=0; r_ind<centered_r_cluster_dmn_t::dmn_size(); r_ind++){
    nfft_R_2_k.f_hat[r_ind][0] = real(F_R(r_ind)); 
    nfft_R_2_k.f_hat[r_ind][1] = imag(F_R(r_ind)); 
  }

  nfft_trafo(&nfft_R_2_k);
  
  for(int k_ind=0; k_ind<target_k_dmn_t::dmn_size(); k_ind++){
      output_ptr[k_ind].real( nfft_R_2_k.f[k_ind][0] ); 
      output_ptr[k_ind].imag( nfft_R_2_k.f[k_ind][1] ); 

//     VECTOR_OPERATIONS::PRINT(target_k_dmn_t::get_elements()[k_ind]);
//     cout << real(output_ptr[k_ind]) << "\t" << imag(output_ptr[k_ind]) << endl;
  }

//   cout << endl;
}

/*!
 *  \brief the centered_r_cluster is in row-major order because of the FFTW, on which NFFT relies!
 */
template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type>
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::initialize_centered_r_cluster()
{
  centered_r_cluster::get_elements().resize(0);
  
  switch(DIMENSION)
    {
    case 1:
      {
	for(int i0=-grid_size[0]/2; i0<grid_size[0]/2; i0++)
	  {
	    std::vector<double> r_vec(DIMENSION, 0);
		
	    r_vec[0] = i0;

	    centered_r_cluster::get_elements().push_back(r_vec);
	  }
      }
      break;

    case 2:
      {
	for(int i1=-grid_size[1]/2; i1<grid_size[1]/2; i1++)
	  {
	    for(int i0=-grid_size[0]/2; i0<grid_size[0]/2; i0++)
	      {
		std::vector<double> r_vec(DIMENSION, 0);
		
		r_vec[0] = i1;
		r_vec[1] = i0;

		centered_r_cluster::get_elements().push_back(r_vec);
	      }
	  }
      }
      break;
  
    case 3:
      {

	for(int i2=-grid_size[2]/2; i2<grid_size[2]/2; i2++)
	  {
	    for(int i1=-grid_size[1]/2; i1<grid_size[1]/2; i1++)
	      {
		for(int i0=-grid_size[0]/2; i0<grid_size[0]/2; i0++)
		  {
		    std::vector<double> r_vec(DIMENSION, 0);
		    
		    r_vec[0] = i2;
		    r_vec[1] = i1;
		    r_vec[2] = i0;
		    
		    centered_r_cluster::get_elements().push_back(r_vec);
		  }
	      }
	  }
      }
      break;

    default:
      throw std::logic_error(__FUNCTION__);
    }

  assert(int(centered_r_cluster::get_elements().size()) == centered_r_cluster::get_size());
}

template<typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_SHAPE S, typename target_dmn_type> 
void wannier_interpolation_kernel<cluster_domain<scalar_type, D, N, MOMENTUM_SPACE, S>, target_dmn_type>::initialize_cut_off()
{
  for(int l=0; l<lies_within_cutoff.size(); l++)
    lies_within_cutoff(l) = 0;

  for(int R_ind=0; R_ind<source_r_dmn_t::dmn_size(); R_ind++){

    std::vector<double> R_vec = source_r_dmn_t::get_elements()[R_ind];

    //std::vector<std::vector<double> > r_vecs = source_r_cluster_type::find_equivalent_vectors_with_minimal_distance_to_origin(R_vec).second;
    std::vector<std::vector<double> > r_vecs = cluster_operations::equivalent_vectors(R_vec, source_r_cluster_type::get_super_basis_vectors());

    for(size_t l=0; l<r_vecs.size(); l++){

      //std::vector<double> r_aff = LDA_r_cluster_type::get_lattice_affine_coordinate(r_vecs[l]);
      //std::vector<double> r_aff = source_r_cluster_type::get_lattice_affine_coordinate(r_vecs[l]);
      std::vector<double> r_aff = VECTOR_OPERATIONS::COORDINATES(r_vecs[l], source_r_cluster_type::get_basis_vectors());//source_r_cluster_type::get_lattice_affine_coordinate(r_vecs[l]);

      for(int R_cen_ind=0; R_cen_ind<centered_r_cluster_dmn_t::dmn_size(); R_cen_ind++)
	if(VECTOR_OPERATIONS::L2_NORM(r_aff, centered_r_cluster_dmn_t::get_elements()[R_cen_ind]) < 1.e-3)
	  lies_within_cutoff(R_cen_ind) += r_vecs.size();
    }
  }

  if(DIMENSION==2){
    for(int l=0; l<lies_within_cutoff.size(); l++){
      if(centered_r_cluster_dmn_t::get_elements()[l][0] == -grid_size[1]/2 
	 or centered_r_cluster_dmn_t::get_elements()[l][1] == -grid_size[0]/2)
	lies_within_cutoff(l) = 0;
    }
  }

  if(DIMENSION==3){
    for(int l=0; l<lies_within_cutoff.size(); l++){
      if(centered_r_cluster_dmn_t::get_elements()[l][0] == -grid_size[2]/2 
	 or centered_r_cluster_dmn_t::get_elements()[l][1] == -grid_size[1]/2
	 or centered_r_cluster_dmn_t::get_elements()[l][2] == -grid_size[0]/2)
	lies_within_cutoff(l) = 0;
    }
  }

//   for(int l=0; l<lies_within_cutoff.size(); l++){
//     if(lies_within_cutoff(l)>0){
//       VECTOR_OPERATIONS::PRINT(centered_r_cluster_dmn_t::get_elements()[l]);
//       cout << "\t" << lies_within_cutoff(l) << endl;
//     }
//   }
//   cout << endl;
}

#endif
