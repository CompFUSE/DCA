// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements a Wannier interpolation technique.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_HPP

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/interpolation/wannier_interpolation/wannier_interpolation_domain_type.hpp"
#include "dca/phys/domains/cluster/interpolation/wannier_interpolation/wannier_interpolation_generic.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename source_dmn_type, typename target_dmn_type>
class wannier_interpolation {
public:
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output) {
    typedef
        typename wannier_interpolation_domain_type<domain_input, source_dmn_type, target_dmn_type>::Result
            wannier_interpolation_domain;

    typedef typename domain_output::this_type domain_output_list_type;
    dca::util::assert_same<domain_output_list_type, wannier_interpolation_domain>();

    typedef typename domain_input::this_type type_list_input;
    typedef typename domain_output::this_type type_list_output;

    // TODO: Move this static_assert to wannier_interpolation_domain_type?
    static_assert(dca::util::IndexOf<source_dmn_type, type_list_input>::value > -1,
                  "Type list error");

    wannier_interpolation_generic<
        type_list_input, type_list_output, source_dmn_type, target_dmn_type, 0,
        dca::util::IndexOf<source_dmn_type, type_list_input>::value>::execute(f_input, f_output);
  }
};

template <typename source_dmn_type, typename target_dmn_type>
class wannier_interpolation<func::dmn_0<source_dmn_type>, func::dmn_0<target_dmn_type>> {
public:
  // Strips dmn_0 off and calls execute of the wannier_interpolation class from above.
  template <typename scalartype_input, class domain_input, typename scalartype_output, class domain_output>
  static void execute(func::function<scalartype_input, domain_input>& f_input,
                      func::function<scalartype_output, domain_output>& f_output) {
    wannier_interpolation<source_dmn_type, target_dmn_type>::execute(f_input, f_output);
  }
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_HPP
