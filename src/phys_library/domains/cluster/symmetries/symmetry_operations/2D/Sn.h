//-*-C++-*-

#ifndef Sn_H_
#define Sn_H_

#include <cmath>
#include "group_action.h"

/*!
 *  \author Peter Staar
 *  \brief  Reflection over an axis which has an angle 2*pi*n/m with the x-axis.
 */
template<int n, int m>
class Sn_2D : public group_action<2>
{
public: 

  typedef group_action<2> base_type;
  typedef Sn_2D<n,m>      this_type;

  Sn_2D(){};
  ~Sn_2D(){};

  static double* matrix()
  {
//     double theta = double(n)/double(m);
//     double c = cos(2.*M_PI*theta);
//     double s = sin(2.*M_PI*theta);
    
//     // in column major order !
//     static double matrix[2*2] = { c*c - s*s, 2*s*c,
// 				  2*s*c    , -(c*c - s*s)};
//     return matrix;

    static double* matrix = init();
    return matrix;
  }

private:

  static double* init()
  {
    static double* matrix = new double[4];

    //double theta = double(n)/double(m);
    double c = COSINE_EVAL<n,m>::value();//cos(2.*M_PI*theta);
    double s = SINE_EVAL  <n,m>::value();//sin(2.*M_PI*theta);

//     cout.precision(16);
//     cout << "\t" << c << "\t" << s <<endl;
    
    // in column major order !
    matrix[0] = c*c-s*s;
    matrix[1] = 2*s*c;
    matrix[2] = 2*s*c;
    matrix[3] = -(c*c-s*s);

    return matrix;
  }

};


#endif
