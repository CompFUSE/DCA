//-*-C++-*-


#ifndef MAKE_G4_MATRIX_H
#define MAKE_G4_MATRIX_H

namespace dca {

  /*! \file MAKE_G4_MATRIX.h
   *
   * author : peter staar
   */
  template<class parameter_type, class MOMS_type>
  class make_G4_matrix 
  {
#include "type_definitions.h"
    
  public:
    
    template<typename scalartype_1, typename scalartype_2>
    static void execute(FUNC_LIB::function<scalartype_1, dmn_8<b,b,b,b,k_DCA,k_DCA,w_VERTEX,w_VERTEX> >&                  G4,
			FUNC_LIB::function<scalartype_2, dmn_2< dmn_4<b,b,k_DCA,w_VERTEX>, dmn_4<b,b,k_DCA,w_VERTEX> > >& G4_matrix)
    {
      int* coor_1 = new int[G4       .signature()];
      int* coor_2 = new int[G4_matrix.signature()];
      
      for(int i=0; i<G4.size(); i++)
	{
	  G4.linind_2_subind(i, coor_2);
	  
	  coor_1[0] = coor_2[0];
	  coor_1[1] = coor_2[1];
	  coor_1[2] = coor_2[4];//k_1
	  coor_1[3] = coor_2[6];//w_1
	  coor_1[4] = coor_2[2];
	  coor_1[5] = coor_2[3];
	  coor_1[6] = coor_2[5];//k_2
	  coor_1[7] = coor_2[7];//w_2
	  
	  G4_matrix(coor_1[0], coor_1[1], coor_1[2], coor_1[3], coor_1[4], coor_1[5], coor_1[6], coor_1[7])
	    = G4(coor_2[0], coor_2[1], coor_2[2], coor_2[3], coor_2[4], coor_2[5], coor_2[6], coor_2[7]);
	}
      
      delete [] coor_1;
      delete [] coor_2;
    }

  };
}

#endif
