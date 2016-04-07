//-*-C++-*-

#ifndef ADVANCED_FERMIONIC_OVERLAP_MATRIX_ELEMENT_H
#define ADVANCED_FERMIONIC_OVERLAP_MATRIX_ELEMENT_H
#include"phys_library/domain_types.hpp"
using namespace types;

namespace DCA
{
  namespace ADVANCED_EXACT_DIAGONALIZATION
  {
    template<typename parameter_type, typename ed_options>//b_dmn, typename s_dmn, typename r_dmn>
    struct sparse_element
    {

//       typedef ED_type_definitions<parameter_type, b_dmn, s_dmn, r_dmn> ED_type_def;

//       typedef typename ED_type_def::scalar_type  scalar_type;
//       typedef typename ED_type_def::complex_type complex_type;

      typedef typename ed_options::scalar_type  scalar_type;
      typedef typename ed_options::complex_type complex_type;

      int i, j;
      complex_type value;
    };

    template<typename parameter_type, typename ed_options>//typename b_dmn, typename s_dmn, typename r_dmn>
    bool operator< (const sparse_element<parameter_type, ed_options/*b_dmn, s_dmn, r_dmn*/>& el_1,
                    const sparse_element<parameter_type, ed_options/*b_dmn, s_dmn, r_dmn*/>& el_2)
    {
      if(el_1.j < el_2.j)
        return true;

      else if(el_1.j == el_2.j && el_1.i < el_2.i)
        return true;

      else
        return false;
    }

  }

}

#endif
