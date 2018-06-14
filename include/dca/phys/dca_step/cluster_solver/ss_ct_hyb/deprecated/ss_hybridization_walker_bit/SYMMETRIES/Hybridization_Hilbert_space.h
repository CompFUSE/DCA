// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.

#ifndef HYBRIDIZATION_HILBERT_SPACE_H
#define HYBRIDIZATION_HILBERT_SPACE_H

namespace QMC 
{
  class Hybridization_Hilbert_space
  {

#include "type_definitions.h"
    
  public:

    static void initialize_Hilbert_space(int symm);

    static void set_quantum_number(int& state_number, symmetry_operation symmetrie_number, int& quantum_nb);

    static void sort();
    static void sort_on(symmetry_operation_type symmetrie);

    static int size();
    static std::vector<int>& sizes();

    static void print();

  private:

    static bool double_occupancy(int* state);

  public:

    static int state_size;
    static int space_size;
    static int symmetries;

    static std::vector<int> order;
    static std::vector<Hybridization_quantum_state> space;
    static std::vector<int> symmetrie_sizes;

  };

  int Hybridization_Hilbert_space::state_size = 0;
  int Hybridization_Hilbert_space::space_size = 0;
  int Hybridization_Hilbert_space::symmetries = 0;

  std::vector<Hybridization_quantum_state>  Hybridization_Hilbert_space::space;
  std::vector<int>                               Hybridization_Hilbert_space::symmetrie_sizes;

  void Hybridization_Hilbert_space::initialize_Hilbert_space(int symm)
  {
    symmetries = symm;
    state_size = b::dmn_size()*r_DCA::dmn_size()*s::dmn_size();
    int max_space_size =  pow(2.,state_size);

    int* state  = new int[state_size];
    int* powers = new int[state_size];

    powers[state_size-1] = 1;
    for(int i = state_size-2; i>-1; i--)
      powers[i] = powers[i+1]*2;

    space_size = 0;
    for(int i = 0; i < max_space_size; i++){
      for(int j = 0; j < state_size; j++){
	state[j] = ((int) i/powers[state_size-1-j])%2;
      }
      //if(!double_occupancy(state)){
	space.push_back(Hybridization_quantum_state(state_size, state, symmetries, i));
	space_size++;
	// }
    }


    delete [] state;
    delete [] powers;
  }

  void Hybridization_Hilbert_space::set_quantum_number(int& state_number, symmetry_operation symmetrie_number, int& quantum_nb)
  {
    space[state_number].set_quantum_number(state_number, symmetrie_number, quantum_nb);
  }

  void Hybridization_Hilbert_space::sort()
  {
    std::sort(space.begin(), space.end());
  }

  void Hybridization_Hilbert_space::sort_on(symmetry_operation_type symmetrie)
  {
    Hybridization_quantum_state::compare_type = symmetrie;
    std::sort(space.begin(), space.end(), Hybridization_quantum_state::sort_on);
  }

  void Hybridization_Hilbert_space::print()
  {
    cout << "\n \t Hilbert space : " << endl;
    for(int i = 0; i < space_size; i++){
      cout << "\t";
      space[i].print();
    }
  }

  bool Hybridization_Hilbert_space::double_occupancy(int* state)
  {
    for(int i = 0; i<state_size/2; i++)
      if(state[i] == 1 && state[i+state_size/2] == 1)
	return true;
    return false;
  }

  int Hybridization_Hilbert_space::size()
  {
    static int size = (int) (sizes().size());
    return size;
  }

  std::vector<int>& Hybridization_Hilbert_space::sizes()
  {
    symmetrie_sizes.resize(1,1);
    for(int i = 1; i< space_size; i++){
      if(space[i].quantum_number == space[i-1].quantum_number)
	 symmetrie_sizes[(int) symmetrie_sizes.size()-1] += 1;
      else
	symmetrie_sizes.push_back(1);
    } 

    return symmetrie_sizes;
  }

}
#endif
