// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Domain base class.

#ifndef DCA_FUNCTION_DOMAINS_DOMAIN_HPP
#define DCA_FUNCTION_DOMAINS_DOMAIN_HPP

#include <vector>

namespace dca {
namespace func {
// dca::func::

class domain {
public:
  domain();

  void reset();

  // domain interface functions
  const int& get_size() const {
    return size;
  }

  const std::vector<int>& get_branch_domain_sizes() {
    return branch_domain_sizes;
  }
  int get_Nb_branch_domains() const {
    return branch_domain_sizes.size();
  }
  int get_branch_size(int branch_index) const {
    return branch_domain_sizes[branch_index];
  }

  const std::vector<int>& get_leaf_domain_sizes() const {
    return leaf_domain_sizes;
  }
  int get_Nb_leaf_domains() const {
    return leaf_domain_sizes.size();
  }
  int get_subdomain_size(int subdomain_index) const {
    return leaf_domain_sizes[subdomain_index];
  }

  const std::vector<int>& get_leaf_domain_steps() const {
    return leaf_domain_steps;
  }
  const std::vector<int>& get_branch_domain_steps() const {
    return branch_domain_steps;
  }

  // linind <--> subdomain_indices != branch_indices
  void linind_2_subind(int linind, int* subind) const;
  void subind_2_linind(const int* subind, int& linind) const;

protected:
  int size;

  std::vector<int> branch_domain_sizes;
  std::vector<int> leaf_domain_sizes;

  std::vector<int> branch_domain_steps;
  std::vector<int> leaf_domain_steps;
};

}  // func
}  // dca

#endif  // DCA_FUNCTION_DOMAINS_DOMAIN_HPP
