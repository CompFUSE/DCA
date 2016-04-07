//-*-C++-*-

/*! \file Hybridization_symmetries.h  
 *
 * author Bart Ydens
 *
 * 
 */

#ifndef HYBRIDIZATION_SYMMETRIES_H
#define HYBRIDIZATION_SYMMETRIES_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace QMC 
{
  template<symmetry_operation_type symmetrie>
  class Hybridization_symmetries
  {


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
