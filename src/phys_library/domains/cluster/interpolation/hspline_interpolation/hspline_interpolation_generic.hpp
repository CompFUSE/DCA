// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the generic loop over all the subdomains through template recursion.

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_INTERPOLATION_HSPLINE_INTERPOLATION_HSPLINE_INTERPOLATION_GENERIC_HPP
#define PHYS_LIBRARY_DOMAINS_CLUSTER_INTERPOLATION_HSPLINE_INTERPOLATION_HSPLINE_INTERPOLATION_GENERIC_HPP

#include "dca/function/function.hpp"
#include "dca/util/type_list.hpp"
#include "phys_library/domains/cluster/interpolation/hspline_interpolation/hspline_interpolation_any_2_any.hpp"

template <typename type_list1, typename type_list2, typename type_input, typename type_output,
          int dmn_shift, int next_index>
struct hspline_interpolation_generic {
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output, double a) {
    // typedef typename TypeListAt<type_list1,IndexOf<type_list1, type_input>::value>::Result
    // new_typelist1;
    // typedef typename TypeListAt<type_list2,IndexOf<type_list1, type_input>::value>::Result
    // new_typelist2;

    hspline_interpolation_any_2_any<type_input, type_output,
                                    dca::util::IndexOf<type_input, type_list1>::value +
                                        dmn_shift>::execute(f_input, f_output, a);
  }
};

// End of recursion: next_index = -1
template <typename type_list1, typename type_list2, typename type_input, typename type_output, int dmn_shift>
struct hspline_interpolation_generic<type_list1, type_list2, type_input, type_output, dmn_shift, -1> {
  template <typename scalartype_1, typename dmn_type_1, typename scalartype_2, typename dmn_type_2>
  static void execute(func::function<scalartype_1, dmn_type_1>& /*f_source*/,
                      func::function<scalartype_2, dmn_type_2>& /*F_target*/, double /*a*/) {}
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_INTERPOLATION_HSPLINE_INTERPOLATION_HSPLINE_INTERPOLATION_GENERIC_HPP
