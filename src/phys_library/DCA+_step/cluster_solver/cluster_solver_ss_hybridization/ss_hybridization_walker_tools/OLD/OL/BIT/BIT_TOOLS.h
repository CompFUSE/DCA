//-*-C++-*-

/*! \file BIT_TOOLS.h  
 *
 *  author Peter Staar 
 */

#ifndef BIT_TOOLS_H
#define BIT_TOOLS_H

namespace QMC {
  
  struct BIT_TOOLS
  {
#include "type_definitions.h"

  public :

    template<typename configuration_t, typename vertex_vertex_matrix_type, typename function_type_2>
    static void FULL_CHECK(configuration_t& configuration, function<vertex_vertex_matrix_type, nu>& M, function_type_2& F);

    template<typename configuration_t, typename hamiltonian_type>
    static double TRACE_CHECK(configuration_t& configuration, hamiltonian_type& Hamiltonian);

  };

  template<typename configuration_t, typename vertex_vertex_matrix_type, typename function_type_2>
  void BIT_TOOLS::FULL_CHECK(configuration_t& configuration, function<vertex_vertex_matrix_type, nu>& M, function_type_2& F)
  {
    typedef vertex_vertex_matrix_type                         v_v_matrix_t;
    typedef typename configuration_t::orbital_configuration_t orbital_configuration_t;

    static int* coor = new int[2];
    static nu nu_obj;

    int N_spin_orbitals = nu_obj.get_size();

    for(int i=0; i<N_spin_orbitals; i++)
      {
	nu_obj.linind_2_subind(i, coor);

	orbital_configuration_t& vertices = configuration.get_vertices(i);

	size_t N_v = vertices.size();

// 	cout << i << "\t <k> = " << N_v << endl;
// 	for(size_t j=0; j<N_v; j++)
// 	  cout << "\t {" << vertices[j].t_start() << " , " << vertices[j].t_end() << " }"; 
// 	cout << "\n";

	assert(M(i).get_current_size() == int(N_v));

	if(N_v>0)
	  {
	    v_v_matrix_t F_new(N_v,N_v);
	    v_v_matrix_t M_new(N_v,N_v);
	    
	    assert(int(N_v) == F_new.get_current_size());
	    assert(int(N_v) == M_new.get_current_size());

	    for(size_t l1=0; l1<N_v; l1++){
	      for(size_t l2=0; l2<N_v; l2++){
		F_new(l1,l2) = HYBRIDIZATION_TOOLS::interpolate_F(coor, vertices[l1].t_end()-vertices[l2].t_start(), F);
		assert(M(i)(l1,l2) == M(i)(l1,l2)); // check for nan's
	      }
	    }

	    invert_plan<double> inv_pln(N_v, F_new.get_global_size());
	    memcpy(inv_pln.Matrix, &F_new(0,0), sizeof(double)*F_new.get_global_size()*F_new.get_global_size());
	    inv_pln.execute_plan();
	    memcpy( &M_new(0,0), inv_pln.inverted_matrix,sizeof(double)*F_new.get_global_size()*F_new.get_global_size());
	    
	    if(!M(i).difference(M_new))
	      {
		cout << "F-matrix :: \n\n";
		F_new.print();
		cout << "reference-matrix :: \n\n";
		M_new.print();
		cout << "real-matrix :: \n\n";
		M(i).print();
		assert(false);
	      }
	  }

// 	cout << "\n";
      }
  }

  template<typename configuration_t, typename hamiltonian_type>
  double BIT_TOOLS::TRACE_CHECK(configuration_t& configuration, hamiltonian_type& Hamiltonian)
  {
    typedef typename configuration_t::orbital_configuration_t      orbital_configuration_t;
    typedef typename hamiltonian_type::lambda_type                 lambda_type;
    typedef typename hamiltonian_type::Hilbert_domain              Hilbert_domain;
    typedef typename hamiltonian_type::creation_annihilation_type  c_a_t;
    typedef typename hamiltonian_type::lambda_type                 lambda_type;

    int N_spin_orbitals = (int) nu::dmn_size();

    std::vector<std::pair<double,int> > full_configuration;

    for(int i=0; i<N_spin_orbitals; i++){
      orbital_configuration_t& vertices = configuration.get_vertices(i);
      for(int j = 0; j<(int)vertices.size(); j++){
	full_configuration.push_back(std::pair<double, int>(vertices[j].t_start(), 2*i));
	full_configuration.push_back(std::pair<double, int>(vertices[j].t_end()  , 2*i+1));
      }
    }

    sort(full_configuration.begin(),full_configuration.end());
    
    double trace = 0;

    lambda_type                      exponentials;
    c_a_t                            matrix1;
    c_a_t                            matrix2;

    block_gemm_plan<std::complex<double> > block_plan;

    for(int i=0; i<N_spin_orbitals; i++){
      orbital_configuration_t& vertices = configuration.get_vertices(i);
      if(vertices.size() > 0){}
      else if(configuration.get_full_line(i)){
	full_configuration.push_back(std::pair<double, int>(time_domain_type::beta, 2*i+1));
	full_configuration.push_back(std::pair<double, int>(time_domain_type::beta, 2*i));
      }
      else{
	full_configuration.push_back(std::pair<double, int>(time_domain_type::beta, 2*i));
	full_configuration.push_back(std::pair<double, int>(time_domain_type::beta, 2*i+1));
      }
    }
    int fc_size = (int)full_configuration.size();

    for(int i = 0; i<Hilbert_domain::dmn_size(); i++){
      exponentials(i,i).resize(Hamiltonian.sizes[i]);
      for(int j = 0; j< Hamiltonian.sizes[i]; j++){
	if(full_configuration[fc_size-1].first != time_domain_type::beta)
	  exponentials(i,i)[j] = std::exp((full_configuration[fc_size-1].first-time_domain_type::beta)*Hamiltonian.lambda(i,i)[j]);
	else
	  exponentials(i,i)[j] = 1;
      }
    }
    //Hamiltonian.print_matrix(Hamiltonian.H);
      
    int type = full_configuration[fc_size-1].second%2;

    if(type == 0)
      block_plan.execute_plan_resize(exponentials, Hamiltonian.c_dagger(full_configuration[fc_size-1].second/2), matrix1);
    else{
      assert(type == 1);
      block_plan.execute_plan_resize(exponentials, Hamiltonian.c_minus (full_configuration[fc_size-1].second/2), matrix1);
    }
    //Hamiltonian.print_matrix(matrix1);
    for(int k = fc_size-1; k>0; k--)
      {
	assert((full_configuration[k].first-full_configuration[k-1].first) >= 0.);

	for(int i = 0; i<Hilbert_domain::dmn_size(); i++)
	  for(int j = 0; j<Hamiltonian.sizes[i]; j++)
	    if(full_configuration[k].first != full_configuration[k-1].first)
	      exponentials(i,i)[j] = std::exp(-(full_configuration[k].first-full_configuration[k-1].first)*Hamiltonian.lambda(i,i)[j]);
	    else
	      exponentials(i,i)[j] = 1;

	block_plan.execute_plan_resize(matrix1,exponentials,matrix2);

	type = full_configuration[k-1].second%2;

	if(type == 0)
	  block_plan.execute_plan_resize(matrix2, Hamiltonian.c_dagger(full_configuration[k-1].second/2), matrix1, 1., 0., 'N', 'N');
	else{
	  assert(type == 1);
	  block_plan.execute_plan_resize(matrix2, Hamiltonian.c_minus (full_configuration[k-1].second/2), matrix1, 1., 0., 'N', 'N');
	}
	//Hamiltonian.print_matrix(matrix1);
      }

    for(int i = 0; i<Hilbert_domain::dmn_size();i++)
      for(int j = 0; j<Hamiltonian.sizes[i];j++)
	exponentials(i,i)[j] = std::exp(-full_configuration[0].first*Hamiltonian.lambda(i,i)[j]);

    block_plan.execute_plan_resize(matrix1,exponentials,matrix2);

    for(int i = 0; i<Hilbert_domain::dmn_size(); i++)
      for(int j = 0; j<Hamiltonian.sizes[i]; j++)
	if(matrix2(i,i).get_current_size() != std::pair<int,int>(0,0)){
	  trace += matrix2(i,i)(j,j).real();
	}
    // Hamiltonian.print_matrix(matrix2);
    cout << "trace : " << trace << endl;
    configuration.print();

    
    return trace;
  }
}

#endif
