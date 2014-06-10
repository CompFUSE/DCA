//-*-C++-*-

#ifndef IDENTITY_GROUP_OPERATION_H
#define IDENTITY_GROUP_OPERATION_H

#include <cmath>
#include "group_action.h"

/*!
 *      Author: Peter Staar
 *
 *  \brief Reflection over an axis which has an angle 2*pi*n/m with the x-axis.
 */

template<int DIMENSION>
class identity_group_operation : public group_action<DIMENSION>
{};

template<>
class identity_group_operation<2> : public group_action<2>
{
public: 

  typedef group_action<2>             base_type;
  typedef identity_group_operation<2> this_type;

  identity_group_operation(){};
  ~identity_group_operation(){};

  const static double* matrix()
  {
    const static double matrix[2*2] = {1., 0.,
				       0., 1.};
    return matrix;
  }
};

template<>
class identity_group_operation<3> : public group_action<3>
{
public: 

  typedef group_action<3>             base_type;
  typedef identity_group_operation<3> this_type;

  identity_group_operation(){};
  ~identity_group_operation(){};

  const static double* matrix()
  {
    const static double matrix[3*3] = {1., 0., 0.,
				       0., 1., 0.,
                                       0., 0., 1.};
    return matrix;
  }
};

#endif
