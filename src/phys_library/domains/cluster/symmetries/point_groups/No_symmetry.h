//-*-C++-*-

#ifndef NULL_SYMMETRY_H
#define NULL_SYMMETRY_H

#include "phys_library/domains/cluster/symmetries/point_group.h"
#include "dca/util/type_list.hpp"
#include "phys_library/domains/cluster/symmetries/symmetry_operations/identity_operation.h"
/*!
 *  \author: Peter Staar
 */
template <int DIMENSION>
struct no_symmetry {
  typedef dca::util::Typelist<identity_group_operation<DIMENSION>> point_group_type_list;
};

#endif
