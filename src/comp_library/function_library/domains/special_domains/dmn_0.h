// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Ultimate leaf domain of the domain tree

#ifndef COMP_LIBRARY_FUNCTION_LIBRARY_DOMAINS_SPECIAL_DOMAINS_DMN_0_H
#define COMP_LIBRARY_FUNCTION_LIBRARY_DOMAINS_SPECIAL_DOMAINS_DMN_0_H

#include "dca/util/type_list.hpp"
#include "comp_library/function_library/domains/domain.h"

template <typename parameters>
class dmn_0 : public domain {
public:
  typedef dca::util::Typelist<parameters> this_type;

  typedef parameters parameter_type;
  typedef typename parameters::element_type element_type;

public:
  dmn_0();
  ~dmn_0();

  void initialize();

  void reset();

  static int dmn_size();

  static std::vector<element_type>& get_elements();

  static void print_2_file(const char* filename);
};

template <typename parameters>
dmn_0<parameters>::dmn_0() : domain() {
  dmn_0::initialize();
}

template <typename parameters>
dmn_0<parameters>::~dmn_0() {}

template <typename parameters>
void dmn_0<parameters>::initialize() {
  size = parameters::get_size();

  branch_domain_sizes.push_back(parameters::get_size());
  leaf_domain_sizes.push_back(parameters::get_size());
}

template <typename parameters>
void dmn_0<parameters>::reset() {
  this->domain::reset();

  size = parameters::get_size();

  branch_domain_sizes.push_back(parameters::get_size());
  leaf_domain_sizes.push_back(parameters::get_size());
}

template <typename parameters>
int dmn_0<parameters>::dmn_size() {
  return parameters::get_size();
}

template <typename parameters>
std::vector<typename parameters::element_type>& dmn_0<parameters>::get_elements() {
  return parameters::get_elements();
}

template <typename parameters>
void dmn_0<parameters>::print_2_file(const char* filename) {
  parameters::print_2_file(filename);
}

#endif  // COMP_LIBRARY_FUNCTION_LIBRARY_DOMAINS_SPECIAL_DOMAINS_DMN_0_H
