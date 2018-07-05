// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Ultimate leaf domain of the domain tree.

#ifndef DCA_FUNCTION_DOMAINS_DMN_0_HPP
#define DCA_FUNCTION_DOMAINS_DMN_0_HPP

#include <string>

#include "dca/function/domains/domain.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace func {
// dca::func::

template <typename parameters>
class dmn_0 : public domain {
public:
  typedef dca::util::Typelist<parameters> this_type;

  typedef parameters parameter_type;
  typedef typename parameters::element_type element_type;

  dmn_0();

  void reset();

  static int dmn_size() {
    return parameters::get_size();
  }

  static std::vector<element_type>& get_elements() {
    return parameters::get_elements();
  }

  static const std::string& get_name() {
    static std::string name = "dmn_0<" + parameters::get_name() + ">";
    return name;
  }

  static void print_2_file(const char* filename) {
    parameters::print_2_file(filename);
  }

private:
  void initialize();
};

template <typename parameters>
dmn_0<parameters>::dmn_0() : domain() {
  dmn_0::initialize();
}

template <typename parameters>
void dmn_0<parameters>::initialize() {
  size = parameters::get_size();

  branch_domain_sizes.push_back(parameters::get_size());
  leaf_domain_sizes.push_back(parameters::get_size());
  branch_domain_steps.push_back(1);
  leaf_domain_steps.push_back(1);
}

template <typename parameters>
void dmn_0<parameters>::reset() {
  this->domain::reset();
  dmn_0::initialize();
}

}  // func
}  // dca

#endif  // DCA_FUNCTION_DOMAINS_DMN_0_HPP
