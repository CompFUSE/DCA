//-*-C++-*-

#ifndef Cn_2D_H_
#define Cn_2D_H_

#include "group_action.h"

/*!
 *  \author Peter Staar
 *
 *  \brief Rotation over a 2*pi*n/m
 */
template<int n, int m>
class Cn_2D : public group_action<2>
{
public:

  typedef group_action<2> base_type;
  typedef Cn_2D<n,m>      this_type;

  Cn_2D(){};
  ~Cn_2D(){};

   static double* matrix()
  {
//     double theta = double(n)/double(m);

//     // in column major order !
//     static double matrix[2*2] = { cos(2.*M_PI*theta), sin(2.*M_PI*theta),
// 				  -sin(2.*M_PI*theta), cos(2.*M_PI*theta)};

    static double* matrix = init();
    return matrix;
  }

private:

  static double* init()
  {
    static double* matrix = new double[4];

    double c = COSINE_EVAL<n,m>::value();//cos(2.*M_PI*theta);
    double s = SINE_EVAL  <n,m>::value();//sin(2.*M_PI*theta);

    // in column major order !
    matrix[0] = c;
    matrix[1] = s;
    matrix[2] =-s;
    matrix[3] = c;

    return matrix;
  }

};

#endif
