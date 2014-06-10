//-*-C++-*-

/*
 *      Author: Peter Staar
 */

// elements of the dihedral groups D_n(m)

#ifndef Sn_3D_H_
#define Sn_3D_H_

#include <cmath>
#include "group_action.h"

template<int axis, int n, int m>
class Sn_3D : public group_action<3>
{};

template<int n, int m>
class Sn_3D<0,n,m> : public group_action<3>
{
public:

  typedef group_action<3>  base_type;
  typedef Sn_3D<0,n,m>     this_type;

  Sn_3D(){};
  ~Sn_3D(){};

  const static double* matrix()
  {
     double c( std::cos(2*M_PI*double(n)/double(m)) );
     double s( std::sin(2*M_PI*double(n)/double(m)) );
    
     static double matrix[3*3] = { 1. , 0.        , 0. ,
					0. , c*c - s*s , 2*s*c,
				 	0. , 2*s*c     , -(c*c - s*s)};

    return matrix;
  }
};

template<int n, int m>
class Sn_3D<1,n,m> : public group_action<3>
{
public:

  typedef group_action<3>  base_type;
  typedef Sn_3D<1,n,m>     this_type;

  Sn_3D(){};
  ~Sn_3D(){};

  const static double* matrix()
  {
     double c( std::cos(2*M_PI*double(n)/double(m)) );
     double s( std::sin(2*M_PI*double(n)/double(m)) );
    
     static double matrix[3*3] = { c*c - s*s, 0. , 2*s*c ,
					0.       , 1. , 0.,
				 	2*s*c    , 0. , -(c*c - s*s)};

    return matrix;
  }
};

template<int n, int m>
class Sn_3D<2,n,m> : public group_action<3>
{
public:

  typedef group_action<3>  base_type;
  typedef Sn_3D<2,n,m>     this_type;

  Sn_3D(){};
  ~Sn_3D(){};

  const static double* matrix()
  {
     double c( std::cos(2*M_PI*double(n)/double(m)) );
     double s( std::sin(2*M_PI*double(n)/double(m)) );
    
     static double matrix[3*3] = { c*c - s*s, 2*s*c        , 0.,
					2*s*c    , -(c*c - s*s) , 0.,
					0.       , 0.           , 1.};

    return matrix;
  }
};
#endif
