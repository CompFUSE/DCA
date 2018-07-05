// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Point group product.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUP_PRODUCT_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUP_PRODUCT_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/product_group_action.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename point_group_1, typename point_group_2>
class point_group_product {};

template <typename point_group_1, typename point_group_2>
class point_group_product_left_2_right {};

// product_group_action<Typelist, Typelist>
template <typename head_1, typename head_2, typename tail_1, typename tail_2>
class point_group_product<dca::util::Typelist<head_1, tail_1>, dca::util::Typelist<head_2, tail_2>> {
public:
  typedef typename point_group_product_left_2_right<
      dca::util::Typelist<head_1, tail_1>, dca::util::Typelist<head_2, tail_2>>::Result product_ij;
  typedef typename point_group_product_left_2_right<
      dca::util::Typelist<head_2, tail_2>, dca::util::Typelist<head_1, tail_1>>::Result product_ji;

  typedef typename dca::util::Append<product_ij, product_ji>::type Result;
};

template <typename head_1, typename head_2, typename tail_1, typename tail_2>
class point_group_product_left_2_right<dca::util::Typelist<head_1, tail_1>,
                                       dca::util::Typelist<head_2, tail_2>> {
public:
  typedef
      typename point_group_product_left_2_right<head_1, dca::util::Typelist<head_2, tail_2>>::Result
          product_0j;
  typedef
      typename point_group_product_left_2_right<tail_1, dca::util::Typelist<head_2, tail_2>>::Result
          product_ij;

  typedef typename dca::util::Append<product_0j, product_ij>::type Result;
};

template <typename head_1, typename head_2, typename tail_2>
class point_group_product_left_2_right<dca::util::Typelist<head_2, tail_2>, dca::util::Typelist<head_1>> {
public:
  typedef
      typename point_group_product_left_2_right<dca::util::Typelist<head_2, tail_2>, head_1>::Result Result;
};

template <typename head_1, typename head_2, typename tail_2>
class point_group_product_left_2_right<dca::util::Typelist<head_1>, dca::util::Typelist<head_2, tail_2>> {
public:
  typedef
      typename point_group_product_left_2_right<head_1, dca::util::Typelist<head_2, tail_2>>::Result Result;
};

template <typename head_1, typename head_2>
class point_group_product_left_2_right<dca::util::Typelist<head_1>, dca::util::Typelist<head_2>> {
public:
  typedef product_group_action<head_1, head_2> product_t;
  typedef dca::util::Typelist<product_t> Result;
};

// product_group_action<T, Typelist>
template <typename head, typename head_2, typename tail_2>
class point_group_product_left_2_right<head, dca::util::Typelist<head_2, tail_2>> {
public:
  typedef typename point_group_product_left_2_right<head, tail_2>::Result new_tail;

  typedef typename dca::util::Swap<dca::util::Typelist<head_2, new_tail>, head_2,
                                   product_group_action<head, head_2>>::Result Result;
};

template <typename head, typename last_head>
class point_group_product_left_2_right<head, dca::util::Typelist<last_head>> {
public:
  typedef product_group_action<head, last_head> product_t;
  typedef dca::util::Typelist<product_t> Result;
};

// product_group_action<Typelist, T>
template <typename head, typename head_1, typename tail_1>
class point_group_product_left_2_right<dca::util::Typelist<head_1, tail_1>, head> {
public:
  typedef typename point_group_product_left_2_right<tail_1, head>::Result new_tail;
  typedef typename dca::util::Swap<dca::util::Typelist<head_1, new_tail>, head_1,
                                   product_group_action<head_1, head>>::Result Result;
};

template <typename head, typename last_head>
class point_group_product_left_2_right<dca::util::Typelist<last_head>, head> {
public:
  typedef product_group_action<last_head, head> product_t;
  typedef dca::util::Typelist<product_t> Result;
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUP_PRODUCT_HPP
