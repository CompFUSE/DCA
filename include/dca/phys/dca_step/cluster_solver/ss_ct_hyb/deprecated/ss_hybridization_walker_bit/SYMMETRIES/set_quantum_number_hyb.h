// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.

#ifndef HYBRIDIZATION_SYMMETRIES_ENUMERATED_H
#define HYBRIDIZATION_SYMMETRIES_ENUMERATED_H

namespace QMC 
{

  template<symmetry_operation_type symmetrie>
  struct set_quantum_number_hyb
  {
  public:
    static void execute()
    {}
  };

  template<>
  struct set_quantum_number_hyb<PARTICLE_NUMBER>
  {
#include "type_definitions.h"

  public:

    static void execute()
    {   
      for(int i = 0; i < Hybridization_Hilbert_space::space_size; i++){
	int n  = 0;
	for(int j = 0; j < Hybridization_Hilbert_space::state_size; j++)
	  n += Hybridization_Hilbert_space::space[i].state[j];
	Hybridization_Hilbert_space::set_quantum_number(i,PARTICLE_NUMBER,n);
      }
    }
  };

  template<>
  struct set_quantum_number_hyb<Sz>
  {
#include "type_definitions.h"

  public:

    static void execute()
    {
      b_r_DCA_s b_r_DCA_s_obj;
      int* coor = new int[3];

      for(int i = 0; i < Hybridization_Hilbert_space::space_size; i++){

	int sz = 0;

	for(int flavor = 0; flavor < Hybridization_Hilbert_space::state_size; flavor++){

	  if(Hybridization_Hilbert_space::space[i].state[flavor] == 1){

	    b_r_DCA_s_obj.linind_2_subind(flavor,coor);

	    sz += (coor[2] == 1 ? 1 : -1);
	  }
	}

	Hybridization_Hilbert_space::set_quantum_number(i, Sz, sz);
      }
    }
  };
  
  template<>
  struct set_quantum_number_hyb<TOTAL_MOMENTUM>
  {
#include "type_definitions.h"

  public:

    static void execute()
    {   
      b_r_DCA_s b_r_DCA_s_obj;
      int* coor = new int[3];

      for(int i = 0; i < Hybridization_Hilbert_space::space_size; i++){

	int k = 0;

	for(int flavor = 0; flavor < Hybridization_Hilbert_space::state_size; flavor++){

	  if(Hybridization_Hilbert_space::space[i].state[flavor] == 1){

	    b_r_DCA_s_obj.linind_2_subind(flavor,coor);

	    k = DCA_k_cluster_type::add(k,coor[1]);
	  }
	}

	Hybridization_Hilbert_space::set_quantum_number(i, TOTAL_MOMENTUM, k);
      }
    }
  };
  
}

#endif
