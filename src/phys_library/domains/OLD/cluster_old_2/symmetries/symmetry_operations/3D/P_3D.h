//-*-C++-*-

#ifndef P_3D_H_
#define P_3D_H_

/*!
 *  \author Peter Staar
 *  \brief  rotation around an axis {ux,uy,uz} with angle th = 2*pi*n/m;
 */
class P_3D : public group_action<3>
{
public:

  typedef group_action<3>  base_type;
  typedef P_3D  this_type;

  P_3D(){};
  ~P_3D(){};

  const static double* matrix()
  {
    static double matrix[3*3] = {-1., 0., 0.,
				 0.,-1., 0.,
				 0., 0.,-1.};
				
    return matrix;
  }
};

#endif
