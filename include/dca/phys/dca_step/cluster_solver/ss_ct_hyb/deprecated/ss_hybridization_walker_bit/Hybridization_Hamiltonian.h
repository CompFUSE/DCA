// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.

#ifndef HYBRIDIZATION_HAMILTONIAN_H
#define HYBRIDIZATION_HAMILTONIAN_H

namespace QMC 
{
  template<class parameters_type, class MOMS_type>
  class Hybridization_Hamiltonian
  {
#include "type_definitions.h"
  public:

    typedef typename parameters_type::Concurrency_Type concurrency_type;

    typedef dmn_0<Hybridization_Hilbert_domain> Hilbert_domain;
    
    typedef FUNC_LIB::function<resizeable_rectangular_matrix<std::complex<double> >, dmn_2<Hilbert_domain,Hilbert_domain> > creation_annihilation_type;
    typedef FUNC_LIB::function<                              std::vector <double>  , dmn_2<Hilbert_domain,Hilbert_domain> > lambda_type;

  public:

    Hybridization_Hamiltonian(parameters_type& parameters_ref,
			      MOMS_type&       MOMS_ref,
			      Hybridization_Hilbert_space& Hilbert_space);
    ~Hybridization_Hamiltonian();

    void initialize();

    template<class stream_type>
    void to_JSON(stream_type& ss);

  private:

    void set_function_sizes();
    void reset_matrices();
   
    void initialize_C();
   
    void initialize_U();
    void initialize_mu();

    void particle_number_coordinate(int kolom, int row, int*& p_nb_coor);

    void initialize_density_H();

    void diagonolize_H_C();

    void initialize_sparse();

  public:

    void print_matrix(FUNC_LIB::function<resizeable_rectangular_matrix<std::complex<double> >,dmn_2<Hilbert_domain,Hilbert_domain> >& matrix);
    
  private:

    parameters_type&  parameters;
    MOMS_type&        MOMS;

    int bdmn;
    int rdmn;
    int sdmn;
    int FLAVORS;

    double LAMBDA_MIN;

  public:

    int* sizes;
  
    FUNC_LIB::function<std::complex<double>, dmn_2<nu,nu> > U; 
    FUNC_LIB::function<std::complex<double>, nu           > mu; 

    FUNC_LIB::function<sparse_block_matrix<resizeable_rectangular_matrix<std::complex<double> >, Hilbert_domain >, nu> c_min;
    FUNC_LIB::function<sparse_block_matrix<resizeable_rectangular_matrix<std::complex<double> >, Hilbert_domain >, nu> c_dag;

    FUNC_LIB::function<FUNC_LIB::function<resizeable_rectangular_matrix<std::complex<double> >, dmn_2<Hilbert_domain,Hilbert_domain> >, nu> c_minus;
    FUNC_LIB::function<FUNC_LIB::function<resizeable_rectangular_matrix<std::complex<double> >, dmn_2<Hilbert_domain,Hilbert_domain> >, nu> c_dagger;
   
    FUNC_LIB::function<resizeable_rectangular_matrix<std::complex<double> >,dmn_2<Hilbert_domain,Hilbert_domain> > H;

    FUNC_LIB::function<resizeable_rectangular_matrix<std::complex<double> >,dmn_2<Hilbert_domain,Hilbert_domain> > V;
    FUNC_LIB::function<                              std::vector <double>  ,dmn_2<Hilbert_domain,Hilbert_domain> > lambda;

    Hybridization_Hilbert_space& Hilbert_space;
  };

  template<class parameters_type, class MOMS_type>
  Hybridization_Hamiltonian<parameters_type, MOMS_type>::Hybridization_Hamiltonian(parameters_type& parameters_ref,
											     MOMS_type&       MOMS_ref,
											     Hybridization_Hilbert_space& Hilb_space):
    parameters(parameters_ref),
    MOMS(MOMS_ref),

    bdmn(b::dmn_size()),
    rdmn(r_DCA::dmn_size()),
    sdmn(s::dmn_size()),
    FLAVORS(bdmn*rdmn*sdmn),

    LAMBDA_MIN(0),

    U("U"),
    mu("mu"),

    c_min("c_min"),
    c_dag("c_dag"),

    c_minus("c_minus"),
    c_dagger("c_dagger"),

    H("HAMILTONIAN"),

    V("V"),
    lambda("lambda"),

    Hilbert_space(Hilb_space)
  {}

  template<class parameters_type, class MOMS_type>
  Hybridization_Hamiltonian<parameters_type, MOMS_type>::~Hybridization_Hamiltonian()
  {
    delete [] sizes;
  }

  template<class parameters_type, class MOMS_type>
  void Hybridization_Hamiltonian<parameters_type, MOMS_type>::initialize()
  {
    set_function_sizes();

    reset_matrices();

    initialize_U();
    initialize_mu();

    initialize_C();

    initialize_density_H();
    diagonolize_H_C();

    initialize_sparse();
    
    // for(int flavor = 0; flavor<FLAVORS; flavor++)
    //c_min(0).print();

    // print_matrix(H);
  }
  
  template<class parameters_type, class MOMS_type>
  void Hybridization_Hamiltonian<parameters_type, MOMS_type>::set_function_sizes()
  {
    sizes = new int[Hilbert_domain::dmn_size()];
    for(int i = 0; i< Hilbert_domain::dmn_size();i++)
      sizes[i] = Hilbert_domain::get_elements()[i];

    for(int k = 0; k<Hilbert_domain::dmn_size();k++)
      {
	H(k,k)     .resize_no_copy(std::pair<int,int>(sizes[k],sizes[k]),
				   std::pair<int,int>(sizes[k],sizes[k]));

	V(k,k)     .resize_no_copy(std::pair<int,int>(sizes[k],sizes[k]),
				   std::pair<int,int>(sizes[k],sizes[k]));

	lambda(k,k).resize        (                            sizes[k]);
      }
  }

  template<class parameters_type, class MOMS_type>
  void Hybridization_Hamiltonian<parameters_type, MOMS_type>::reset_matrices()
  {
    U.reset();
    mu.reset();

    c_minus.reset();
    c_dagger.reset();
  }


  template<class parameters_type, class MOMS_type>  
  void Hybridization_Hamiltonian<parameters_type, MOMS_type>::initialize_U()
  {
    for(int i=0; i<nu::dmn_size(); i++){
      for(int j=0; j<nu::dmn_size(); j++){
	U(i,j) = MOMS.H_interactions(i,j,0);
      }
    }
  }

 template<class parameters_type, class MOMS_type>  
  void Hybridization_Hamiltonian<parameters_type, MOMS_type>::initialize_mu()
  {
    for(int i=0; i<nu::dmn_size(); i++){
      mu(i) = parameters.get_chemical_potential();
      for(int j=0; j<nu::dmn_size(); j++)
      	mu(i) += (1./2.)*(MOMS.H_interactions(j,i,0)+MOMS.H_interactions(i,j,0))*(1./2.);
    }
  }

  template<class parameters_type, class MOMS_type>  
  void Hybridization_Hamiltonian<parameters_type, MOMS_type>::initialize_C()
  { 
    int* p_nb_coor   = new int[4];
    int* powers      = new int[FLAVORS];
    for(int i = 0; i<FLAVORS; i++)
      powers[i] = Hilbert_space.space[i+1].permutation_number;

    int* permutation = new int[int(pow(2.,FLAVORS))];
    for(int i = 0 ; i<Hybridization_Hilbert_space::space_size; i++)
      permutation[Hilbert_space.space[i].permutation_number] = i;

    for(int flavor = 0; flavor < FLAVORS; flavor++){
      for(int i = 0; i<Hybridization_Hilbert_space::space_size; i++){
	if(Hilbert_space.space[i].state[flavor] == 1){
	  particle_number_coordinate(permutation[Hilbert_space.space[i].permutation_number-powers[flavor]],i,p_nb_coor);

	  if(c_minus (flavor)(p_nb_coor[0],p_nb_coor[1]).get_current_size() == std::pair<int,int>(0,0)){
	    c_minus (flavor)(p_nb_coor[0],p_nb_coor[1]).resize_no_copy(std::pair<int,int>(sizes[p_nb_coor[0] ], sizes[p_nb_coor[1] ]),
								       std::pair<int,int>(sizes[p_nb_coor[0] ], sizes[p_nb_coor[1] ]));
	    c_dagger(flavor)(p_nb_coor[1],p_nb_coor[0]).resize_no_copy(std::pair<int,int>(sizes[p_nb_coor[1] ], sizes[p_nb_coor[0] ]),
								       std::pair<int,int>(sizes[p_nb_coor[1] ], sizes[p_nb_coor[0] ]));
	   }
	   c_minus (flavor)(p_nb_coor[0],p_nb_coor[1])(p_nb_coor[2],p_nb_coor[3]) = 1;
	   c_dagger(flavor)(p_nb_coor[1],p_nb_coor[0])(p_nb_coor[3],p_nb_coor[2]) = 1;
	 }
       }
     }

     delete [] permutation;
     delete [] p_nb_coor;
   }

  template<class parameters_type, class MOMS_type>  
  void Hybridization_Hamiltonian<parameters_type, MOMS_type>::initialize_density_H()
  {
    FUNC_LIB::function<resizeable_rectangular_matrix<std::complex<double> >, dmn_2<Hilbert_domain,Hilbert_domain> > tmp1;
    FUNC_LIB::function<resizeable_rectangular_matrix<std::complex<double> >, dmn_2<Hilbert_domain,Hilbert_domain> > tmp2;
    block_gemm_plan<std::complex<double> > block;

    // --> CHEMICAL POTENTIAL
    for(int i = 0; i<nu::dmn_size(); i++)
      if(abs(mu(i))>1e-6)
	for(int s = 0; s<sdmn; s++)
	  block.execute_plan(c_dagger(i), c_minus(i), H, -mu(i)*1./2., 1);

    // --> INTERACTION TERM
    for(int i = 0; i<nu::dmn_size(); i++)
      for(int j = 0; j<nu::dmn_size(); j++)
	if(abs(U(i,j)) > 1e-6){
	  block.execute_plan_resize(c_dagger(i), c_minus(i), tmp1,1., 0., 'N', 'N');
	  block.execute_plan_resize(c_dagger(j), c_minus(j), tmp2,1., 0., 'N', 'N');
	  block.execute_plan(tmp1, tmp2, H, U(i,j)*1./2., 1);
	}
  }

  template<class parameters_type, class MOMS_type>  
  void Hybridization_Hamiltonian<parameters_type, MOMS_type>::diagonolize_H_C()
  {
    FUNC_LIB::function<resizeable_rectangular_matrix<std::complex<double> >, dmn_2<Hilbert_domain, Hilbert_domain> > TMP;

    block_eigensystem_plan<std::complex<double>, HERMITIAN> block_eigen(sizes);
    block_gemm_plan       <std::complex<double>           > block;

    block_eigen.execute_plan(H,V,lambda);

    LAMBDA_MIN = lambda(0,0)[0];
    for(int i = 0; i < Hilbert_domain::dmn_size(); i++)
      for(int j = 0; j<sizes[i]; j++)
    	LAMBDA_MIN = (LAMBDA_MIN>lambda(i,i)[j])? lambda(i,i)[j] : LAMBDA_MIN;

    for(int i = 0; i < Hilbert_domain::dmn_size(); i++){
      for(int j = 0; j<sizes[i]; j++){
    	lambda(i,i)[j] -= LAMBDA_MIN;
    	assert(lambda(i,i)[j] > -1.e-6);
      }
    }


    for(int flavor = 0; flavor < FLAVORS; flavor++)
      {	
    	block.execute_plan_resize(V  , c_dagger(flavor), TMP             , 1., 0., 'C', 'N');
    	block.execute_plan_resize(TMP, V               , c_dagger(flavor), 1., 0., 'N', 'N');
	
    	block.execute_plan_resize(V  , c_minus (flavor), TMP            , 1., 0., 'C', 'N');
    	block.execute_plan_resize(TMP, V               , c_minus(flavor), 1., 0., 'N', 'N');
      }

  }


  template<class parameters_type, class MOMS_type>
  void Hybridization_Hamiltonian<parameters_type, MOMS_type>::initialize_sparse()
  {
    for(int k = 0; k<FLAVORS;k++){
      c_min(k).blocks.clear();
      c_dag(k).blocks.clear();
      for(int i = 0; i<Hilbert_domain::dmn_size(); i++){
	for(int j = 0; j<Hilbert_domain::dmn_size(); j++){
	  if(c_minus(k)(i,j).get_current_size() != std::pair<int,int>(0,0))
	    c_min(k).add_new_empty(std::pair<int,int>(i,j), std::pair<int,int>(sizes[i],sizes[j]));	
	  if(c_dagger(k)(i,j).get_current_size() != std::pair<int,int>(0,0))
	    c_dag(k).add_new_empty(std::pair<int,int>(i,j), std::pair<int,int>(sizes[i],sizes[j]));	
	}
      }
      c_min(k).copy_from(c_minus(k));
      c_dag(k).copy_from(c_dagger(k));	
    }


  }

  template<class parameters_type, class MOMS_type>  
  void Hybridization_Hamiltonian<parameters_type, MOMS_type>::print_matrix(FUNC_LIB::function<resizeable_rectangular_matrix<std::complex<double> >,dmn_2<Hilbert_domain,Hilbert_domain> >& matrix)
  {

    // for(int i = 0; i<Hilbert_domain::dmn_size(); i++)
    //   for(int j = 0; j<Hilbert_domain::dmn_size(); j++)
    //  	if(matrix(i,j).get_current_size() != std::pair<int,int>(0,0)){
    //   	  cout << i << " " << j<< endl;
    // 	  cout << endl;
    // 	  matrix(i,j).print();
    // 	}

    cout << "\n " << matrix.get_name() << endl;
    for(int i = 0; i<Hilbert_domain::dmn_size(); i++){
      for(int l1=0; l1<sizes[i]; l1++){
    	for(int j = 0; j<Hilbert_domain::dmn_size(); j++){
    	  for(int l2=0; l2<sizes[j]; l2++)
    	    {
    	    if(matrix(i,j).get_current_size() != std::pair<int,int>(0,0))
    	      {
    		assert(matrix(i,j).get_current_size().first == sizes[i]);
    		assert(matrix(i,j).get_current_size().second == sizes[j]);
    		cout << "\t" << matrix(i,j)(l1,l2).real();
    	      }
    	    else
    	      cout << "\t" << 0;
    	    }
    	}
    	cout << endl;
      }
    }
  }

  template<class parameters_type, class MOMS_type>
  void Hybridization_Hamiltonian<parameters_type, MOMS_type>::particle_number_coordinate(int kolom, int row, int*& coor)
  {
    int nb = 0;
    int k = kolom;
    int r = row;

    for(int i = 0;i<4;i++)
      coor[i] = 0;

    while(k - sizes[nb] >= 0){
      coor[0]++;
      k -= sizes[nb];
      nb++;
    }

    nb = 0;
    while(r - sizes[nb] >= 0){
      coor[1]++;
      r -= sizes[nb];
      nb++;
    }

    coor[2] = k;
    coor[3] = r;
  }
}

#endif





















