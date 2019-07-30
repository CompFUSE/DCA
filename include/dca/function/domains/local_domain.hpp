// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements a class to store a subset of a domain distributed among different MPI ranks.
//
// INTERNAL: This class can only be used with dmn_variadic when it's wrapped with dmn_0.

#ifndef DCA_FUNCTION_DOMAINS_LOCAL_DOMAIN_HPP
#define DCA_FUNCTION_DOMAINS_LOCAL_DOMAIN_HPP

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "dca/util/integer_division.hpp"

namespace dca {
namespace func {
// dca::func::

template <typename BaseDomain, int id = 0>
class LocalDomain {
public:
  using element_type = typename BaseDomain::element_type;  // For putting LocalDomain in dmn_0.
  using this_type = LocalDomain<BaseDomain>;

  template <class Grouping>
  static void initialize(const Grouping& group);

  static std::size_t get_physical_size() {
    return elements_.size();
  }
  static std::size_t get_offset() {
    return offset_;
  }

  static std::size_t get_size() {
    assert(initialized_);
    return padded_size_;
  }
  static std::string get_name() {
    return "Local_" + BaseDomain::get_name();
  }
  static const std::vector<element_type>& get_elements() {
    assert(initialized_);
    return elements_;
  }
  static bool is_initialized() {
    return initialized_;
  }

private:
  static bool initialized_;
  static std::vector<element_type> elements_;
  static std::size_t padded_size_;
  static std::size_t offset_;
};

template <class BaseDomain, int id>
bool LocalDomain<BaseDomain, id>::initialized_ = false;
template <class BaseDomain, int id>
std::vector<typename LocalDomain<BaseDomain, id>::element_type> LocalDomain<BaseDomain, id>::elements_;
template <class BaseDomain, int id>
std::size_t LocalDomain<BaseDomain, id>::padded_size_ = 0;
template <class BaseDomain, int id>
std::size_t LocalDomain<BaseDomain, id>::offset_ = 0;

template <class BaseDomain, int id>
template <class Grouping>
void LocalDomain<BaseDomain, id>::initialize(const Grouping& group) {
  if (initialized_) {
    throw(std::logic_error("Domain " + get_name() + " is already initialized."));
  }

  const std::size_t global_size = BaseDomain::get_size();
  padded_size_ = dca::util::ceilDiv(global_size, std::size_t(group.get_size()));

  const std::size_t start = std::min(padded_size_ * group.get_id(), global_size);
  const std::size_t end = std::min(start + padded_size_, global_size);

  const auto physical_size = end - start;
  offset_ = start;
  elements_.resize(physical_size);
  std::copy_n(BaseDomain::get_elements().data() + start, physical_size, elements_.data());

  initialized_ = true;
}

}  // namespace func
}  // namespace dca

#endif  // DCA_FUNCTION_DOMAINS_LOCAL_DOMAIN_HPP
