// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.

#ifndef HYBRIDIZATION_QUANTUM_STATE_H
#define HYBRIDIZATION_QUANTUM_STATE_H

namespace QMC 
{
  class Hybridization_quantum_state
  {
    
  public:

    Hybridization_quantum_state(int state_size, int* state, int symmetries, int permutation_order);
    
    ~Hybridization_quantum_state();

    void set_quantum_number(int& state_number, symmetry_operation symmetrie_number, int& quantum_nb);

    void print();

    static bool sort_on(const Hybridization_quantum_state& state1, const Hybridization_quantum_state& state2);

    bool double_occupancy();

  public:

    int state_size;
    int symmetries;

    std::vector<int> state;
    std::vector<int> quantum_number;

    int permutation_number;

  static symmetry_operation_type compare_type;
  };
  
  symmetry_operation_type Hybridization_quantum_state::compare_type = PARTICLE_NUMBER;
  
  Hybridization_quantum_state::Hybridization_quantum_state(int state_s, int* st, int sym, int perm_number):
    state_size(state_s),
    symmetries(sym),

    state(state_size),
    quantum_number(sym),
    permutation_number(perm_number)
  {
    for(int i =0; i < state_size; i++)
      state[i] = st[i];
  }

  Hybridization_quantum_state::~Hybridization_quantum_state()
  {}

  void Hybridization_quantum_state::set_quantum_number(int& state_number, symmetry_operation symmetrie_number, int& quantum_nb)
  {
    quantum_number[(int) symmetrie_number] = quantum_nb;
  }

  void Hybridization_quantum_state::print()
  {
    cout << "{ { " << state[0];
    for(int i = 1; i < state_size; i++)
      cout << ", "<< state[i] ;
    cout << "}, {" << quantum_number[0];
    for(int i = 1; i < symmetries; i++)
      cout  << ", " << quantum_number[i];
    cout << "}, " << permutation_number << " }" << endl;
  }

  bool Hybridization_quantum_state::sort_on(const Hybridization_quantum_state& state1, const Hybridization_quantum_state& state2)
  {
    if(state1.quantum_number[(int) compare_type] < state2.quantum_number[(int) compare_type])
      return true;
    
    return false;
  }

  bool Hybridization_quantum_state::double_occupancy()
  {
    for(int i = 0; i<state_size;i++)
      if(state[i] == 1 && state[i+state_size/2] == 1)
	return true;

    return false;
  }
  bool operator<(Hybridization_quantum_state state1, Hybridization_quantum_state state2)
  {
    for(int i = 0; i < state1.symmetries; i++){
      if(state1.quantum_number[i] < state2.quantum_number[i])
	return true;
      if(state1.quantum_number[i] > state2.quantum_number[i])
	return false;
    }
    if(state1.permutation_number < state2.permutation_number)
      return true;
    return false;
  }

}
#endif
