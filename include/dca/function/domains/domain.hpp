// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
  std::size_t get_size() const {
    return size;
  }

  const std::vector<std::size_t>& get_branch_domain_sizes() const {
    return branch_domain_sizes;
  }
  int get_Nb_branch_domains() const {
    return branch_domain_sizes.size();
  }
  int get_branch_size(const int branch_index) const {
    return branch_domain_sizes[branch_index];
  }

  const std::vector<std::size_t>& get_leaf_domain_sizes() const {
    return leaf_domain_sizes;
  }
  std::size_t get_Nb_leaf_domains() const {
    return leaf_domain_sizes.size();
  }
  std::size_t get_subdomain_size(const int subdomain_index) const {
    return leaf_domain_sizes[subdomain_index];
  }

  const std::vector<std::size_t>& get_leaf_domain_steps() const {
    return leaf_domain_steps;
  }
  const std::vector<std::size_t>& get_branch_domain_steps() const {
    return branch_domain_steps;
  }

  // linind <--> subdomain_indices != branch_indices
  void linind_2_subind(std::size_t linind, int* subind) const;
  void subind_2_linind(const int* subind, std::size_t& linind) const;

protected:
  std::size_t size;

  std::vector<std::size_t> branch_domain_sizes;
  std::vector<std::size_t> leaf_domain_sizes;

  std::vector<std::size_t> branch_domain_steps;
  std::vector<std::size_t> leaf_domain_steps;
};

}  // namespace func
}  // namespace dca

#endif  // DCA_FUNCTION_DOMAINS_DOMAIN_HPP
