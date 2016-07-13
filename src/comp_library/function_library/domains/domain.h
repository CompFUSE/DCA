// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef COMP_LIBRARY_FUNCTION_LIBRARY_DOMAINS_DOMAIN_H
#define COMP_LIBRARY_FUNCTION_LIBRARY_DOMAINS_DOMAIN_H

#include <cassert>
#include <vector>

class domain {
public:
  domain();

  void reset();

public:
  // domain interface functions
  int& get_size();

  std::vector<int>& get_branch_domain_sizes();
  int get_Nb_branch_domains();
  int get_branch_size(int brach_index);

  std::vector<int>& get_leaf_domain_sizes();
  int get_Nb_leaf_domains();
  int get_subdomain_size(int subdomain_index);

  std::vector<int>& get_leaf_domain_steps();
  std::vector<int>& get_branch_domain_steps();

public:
  // linind <--> subdomain_indices != branch_indices
  void linind_2_subind(int linind, int* subind);
  void subind_2_linind(int* subind, int& linind);

protected:
  int size;

  std::vector<int> branch_domain_sizes;
  std::vector<int> leaf_domain_sizes;

  std::vector<int> branch_domain_steps;
  std::vector<int> leaf_domain_steps;
};

//++++++++++++++++++++++++++++++++++++++++//
//+++  CONSTRUCTOR & DESTRUCTOR        +++//
//++++++++++++++++++++++++++++++++++++++++//

domain::domain()
    : size(0),
      branch_domain_sizes(0),
      leaf_domain_sizes(0),

      branch_domain_steps(0),
      leaf_domain_steps(0) {}

void domain::reset() {
  size = 0;

  leaf_domain_sizes.resize(0);
  branch_domain_sizes.resize(0);

  leaf_domain_steps.resize(0);
  branch_domain_steps.resize(0);
}

//++++++++++++++++++++++++++++++++++++++++//
//+++  DATA-INTERFACE                  +++//
//++++++++++++++++++++++++++++++++++++++++//

int& domain::get_size() {
  return size;
}

std::vector<int>& domain::get_branch_domain_sizes() {
  return branch_domain_sizes;
}

int domain::get_Nb_branch_domains() {
  return branch_domain_sizes.size();
}

int domain::get_branch_size(int branch_index) {
  return branch_domain_sizes[branch_index];
}

std::vector<int>& domain::get_leaf_domain_sizes() {
  return leaf_domain_sizes;
}

int domain::get_Nb_leaf_domains() {
  return leaf_domain_sizes.size();
}

int domain::get_subdomain_size(int subdomain_index) {
  return leaf_domain_sizes[subdomain_index];
}

std::vector<int>& domain::get_leaf_domain_steps() {
  return leaf_domain_steps;
}

std::vector<int>& domain::get_branch_domain_steps() {
  return branch_domain_steps;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++//
//+++  LINEAR INDEX <=> SUB INDEX (subdomains)  +++//
//+++++++++++++++++++++++++++++++++++++++++++++++++//

void domain::linind_2_subind(int linind, int* subind) {
  assert(linind >= 0 && linind < size);

  for (size_t i = 0; i < leaf_domain_sizes.size(); i++) {
    subind[i] = linind % leaf_domain_sizes[i];
    linind = (linind - subind[i]) / leaf_domain_sizes[i];
  }

  assert(linind == 0);
}

void domain::subind_2_linind(int* subind, int& linind) {
  linind = 0;

  for (int i = leaf_domain_sizes.size() - 1; i >= 0; i--) {
    assert(subind[i] < leaf_domain_sizes[i]);
    linind = subind[i] + linind * leaf_domain_sizes[i];
  }

  assert(linind >= 0 && linind < size);
}

#endif  // COMP_LIBRARY_FUNCTION_LIBRARY_DOMAINS_DOMAIN_H
