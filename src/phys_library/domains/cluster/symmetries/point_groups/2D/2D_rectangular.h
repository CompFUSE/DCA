//-*-C++-*-

/*
 *      Author: Peter Staar
 */


#ifndef RECTANGULAR_2D_H_
#define RECTANGULAR_2D_H_

/*
 *  group-actions
 */

  typedef Sn_2D<0,4> Sn_2D_0_4_type;
  typedef Sn_2D<1,4> Sn_2D_1_4_type;

/*
 *  pointgroup :: set of group-actions
 */

typedef Typelist<Sn_2D_0_4_type> Sx_2d;
typedef Typelist<Sn_2D_1_4_type> Sy_2d;

typedef Append<C2, Append<Sx_2d, Sy_2d>::type >::type D2;


#endif
