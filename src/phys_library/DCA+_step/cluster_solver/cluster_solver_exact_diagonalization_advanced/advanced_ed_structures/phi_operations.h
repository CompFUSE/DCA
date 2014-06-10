//-*-C++-*-

#ifndef ADVANCED_FERMIONIC_PHI_OPERATIONS_H
#define ADVANCED_FERMIONIC_PHI_OPERATIONS_H

namespace DCA
{
  namespace ADVANCED_EXACT_DIAGONALIZATION
  {

    template<typename parameter_type, typename ed_options, phi_names phi_name>
    bool operator< (const phi_state<parameter_type, ed_options, phi_name>& phi_obj1, 
		    const phi_state<parameter_type, ed_options, phi_name>& phi_obj2)
    {
      return phi_obj1.phi.to_ulong() < phi_obj2.phi.to_ulong();
    }

    template<typename parameter_type, typename ed_options>
    bool operator== (const phi_state<parameter_type, ed_options, PHI_SINGLET>& phi_obj1, 
		     const phi_state<parameter_type, ed_options, PHI_SINGLET>& phi_obj2)
    {
      return (phi_obj1.phi == phi_obj2.phi and abs(phi_obj1.alpha-phi_obj2.alpha)<ed_options::get_epsilon());
    }
    
  }
}

#endif
