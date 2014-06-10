//-*-C++-*-

#ifndef CLUSTER_REDUCTION_H
#define CLUSTER_REDUCTION_H

/*!
 *  \defgroup SYMMETRIES
 *  \ingroup  ALGORITHMS
 */

// #include "apply_symmetry.h"
// #include "apply_symmetries.h"


/*!
 *  \ingroup SYMMETRIES
 *
 *  \author  Peter Staar
 */
template<class base_cluster_type, class point_group>
class cluster_reduction
{

public:

  cluster_reduction();
  ~cluster_reduction();

  void execute();

private:
};

template<class base_cluster_type, class point_group>
cluster_reduction<base_cluster_type, point_group>::cluster_reduction()
{}

template<class base_cluster_type, class point_group>
cluster_reduction<base_cluster_type, point_group>::~cluster_reduction()
{}

template<class base_cluster_type, class point_group>
void cluster_reduction<base_cluster_type, point_group>::execute()
{
  //cout << __FUNCTION__ << endl;

  search_maximal_symmetry_group<base_cluster_type, point_group, UNIT_CELL>::execute();

  search_maximal_symmetry_group<base_cluster_type, point_group, SUPER_CELL>::execute();

  set_symmetry_matrices<base_cluster_type>::execute();

//   set_symmetry_matrices<base_cluster_type>::print_on_shell();
//   throw std::logic_error(__FUNCTION__);
}

#endif
