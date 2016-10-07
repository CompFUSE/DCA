// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements a Hermite spline interpolation technique.

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_INTERPOLATION_HSPLINE_INTERPOLATION_HSPLINE_INTERPOLATION_HPP
#define PHYS_LIBRARY_DOMAINS_CLUSTER_INTERPOLATION_HSPLINE_INTERPOLATION_HSPLINE_INTERPOLATION_HPP

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/util/type_list.hpp"
#include "phys_library/domains/cluster/interpolation/hspline_interpolation/hspline_interpolation_domain_type.hpp"
#include "phys_library/domains/cluster/interpolation/hspline_interpolation/hspline_interpolation_generic.hpp"

template <typename source_dmn_type, typename target_dmn_type>
class hspline_interpolation {
public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output, double a) {
    typedef
        typename hspline_interpolation_domain_type<domain_input, source_dmn_type, target_dmn_type>::Result
            hspline_interpolation_domain;

    typedef typename domain_output::this_type domain_output_list_type;
    dca::util::assert_same<domain_output_list_type, hspline_interpolation_domain>();

    typedef typename domain_input::this_type type_list_input;
    typedef typename domain_output::this_type type_list_output;

    static_assert(dca::util::IndexOf<source_dmn_type, type_list_input>::value > -1,
                  "Type list error");

    hspline_interpolation_generic<
        type_list_input, type_list_output, source_dmn_type, target_dmn_type, 0,
        dca::util::IndexOf<source_dmn_type, type_list_input>::value>::execute(f_input, f_output, a);
  }
};

template <typename source_dmn_type, typename target_dmn_type>
class hspline_interpolation<func::dmn_0<source_dmn_type>, func::dmn_0<target_dmn_type>> {
public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output, double a) {
    hspline_interpolation<source_dmn_type, target_dmn_type>::execute(f_input, f_output, a);
  }
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_INTERPOLATION_HSPLINE_INTERPOLATION_HSPLINE_INTERPOLATION_HPP
