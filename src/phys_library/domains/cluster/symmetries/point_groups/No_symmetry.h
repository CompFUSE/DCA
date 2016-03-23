//-*-C++-*-

#ifndef NULL_SYMMETRY_H
#define NULL_SYMMETRY_H

#include "point_group.h"

/*!
 *  \author: Peter Staar
 */
template<int DIMENSION>
struct no_symmetry
{
  typedef Typelist<identity_group_operation<DIMENSION>> point_group_type_list;
};

#endif
