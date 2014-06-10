//-*-C++-*-

#ifndef ADVANCED_FERMIONIC_RELATIONAL_OPERATIONAL_H
#define ADVANCED_FERMIONIC_RELATIONAL_OPERATIONAL_H

namespace DCA
{
  namespace ADVANCED_EXACT_DIAGONALIZATION
  {

    template<typename parameter_type, typename ed_options>
      bool operator< (const rep_struct<parameter_type, ed_options>& phi_obj1,
                      const rep_struct<parameter_type, ed_options>& phi_obj2)
      {
        return (phi_obj1.phi.to_ulong() < phi_obj2.phi.to_ulong());
      }

  }

}

#endif
