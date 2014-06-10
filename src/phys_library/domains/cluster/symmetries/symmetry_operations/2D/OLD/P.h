//-*-C++-*-

/*
 *      Author: peterstaar
 */

// elements of the dihedral groups D_n(m)
/*
#ifndef P_2D_H_
#define P_2D_H_

#include <cmath>
#include "group_action.h"


class P_2D : public group_action<2>
{
public:

  typedef group_action<2> base_type;
  typedef P_2D this_type;

  P_2D(){};
  ~P_2D(){};

  const static double* matrix()
  {
    const static double matrix[2*2] = { -1.,  0.,
				         0., -1.};
    return matrix;
  }
};


#endif
*/
