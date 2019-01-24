// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.

#ifndef HYBRIDIZATION_SYMMETRIES_H
#define HYBRIDIZATION_SYMMETRIES_H

namespace QMC 
{
  template<symmetry_operation_type symmetrie>
  class Hybridization_symmetries
  {

#include "type_definitions.h"

  public:

    static void execute_symmetrie();

  private:

    static void set_order(Hybridization_Hilbert_space& Hilbert_space);

    static int least_commen_multiple(int a,int b);
  };
  
  template<symmetry_operation_type symmetrie>
  void Hybridization_symmetries<symmetrie>::execute_symmetrie()
  {
    set_quantum_number_hyb<symmetrie>::execute();
  }

  template<symmetry_operation_type symmetrie>
  int Hybridization_symmetries<symmetrie>::least_commen_multiple(int a,int b)
  {
    int n;
    for(n=1;;n++)
      {
  	if(n%a == 0 && n%b == 0)
  	  return n;
      }
  }

}

#endif























 // template<symmetry_operation_type symmetrie>
  // void Hybridization_full_symmetries<symmetrie>::set_order(Hybridization_full_Hilbert_space& Hilbert_space)
  // {
  //   int order = 1;
  //   for(int i = 0; i < Hilbert_space.space_size; i++){
  //     int j = i;
  //     int current_order = 1;
  //     while(i != Hilbert_space.space[j].permutation_number){
  // 	j = Hilbert_space.space[j].permutation_number;
  // 	current_order++;
  //     }
  //     order = least_commen_multiple(order, current_order);
  //   }
  //   Hilbert_space.order[(int) symmetrie] = order;
  // }
