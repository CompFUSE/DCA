//-*-C++-*-

#ifndef FERMIONIC_STATE_H
#define FERMIONIC_STATE_H

/*!
 *   \author Peter Staar
 */
/*
template<typename b_dmn, typename s_dmn, typename r_dmn>
class fermionic_state
{
public:

  typedef fermionic_operator<b_dmn, s_dmn, r_dmn> f_operator;
  typedef fermionic_state   <b_dmn, s_dmn, r_dmn> f_state;

public:
  
  fermionic_state();

  ~fermionic_state();
  
  int get_N();
  int get_Sz();

public:

  std::vector<f_operator> state;
};

template<typename b_dmn, typename s_dmn, typename r_dmn>
fermionic_state<b_dmn, s_dmn, r_dmn>::fermionic_state():
  state(0)
{}

template<typename b_dmn, typename s_dmn, typename r_dmn>
fermionic_state<b_dmn, s_dmn, r_dmn>::~fermionic_state()
{}

template<typename b_dmn, typename s_dmn, typename r_dmn>
int fermionic_state<b_dmn, s_dmn, r_dmn>::get_N()
{
  return state.size();
}

template<typename b_dmn, typename s_dmn, typename r_dmn>
int fermionic_state<b_dmn, s_dmn, r_dmn>::get_Sz()
{
  int Sz=0;
  for(size_t l=0; l<state.size(); l++)
    Sz += state[l].Sz;

  return Sz; 
}
*/

#endif
