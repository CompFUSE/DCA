// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements a class to store a subset of a previously defined one-dimensional domain.
// This subset can either be a set of single elements or a range starting from either the beginning
// or the middle index (more precisely base-domain-size/2) of the base domain. The domain is created
// by a type definition followed by a call to an initialize function.
//
// INTERNAL: This class is probably good enough for its purpose, but it doesn't seem to be
// 'complete'.

#ifndef COMP_LIBRARY_FUNCTION_LIBRARY_DOMAINS_SPECIAL_DOMAINS_REDUCED_DOMAIN_HPP
#define COMP_LIBRARY_FUNCTION_LIBRARY_DOMAINS_SPECIAL_DOMAINS_REDUCED_DOMAIN_HPP

#include <cassert>
#include <cstdlib>  // for std::size_t
#include <stdexcept>
#include <string>
#include <vector>
#include "comp_library/function_library/domains/special_domains/dmn_variadic.h"

namespace dca {
namespace func {
// dca::func::

template <typename BaseDomain>
class ReducedDomain {
public:
  using ElementType = typename BaseDomain::element_type;

  // If from_middle = true, creates a domain that contains the first n elements of the base domain
  // counted from the middle index (= base-domain-size/2) of the base domain.
  // If from_middle = false, creates a domain that contains the first n elements of the base domain.
  static void initialize(const std::size_t n, const bool from_middle, const std::string& name = "");

  // Creates a domain that only contains the elements of the base domain with the given indices.
  // INTERNAL: What should be the preconditions on indicies? Sorted, no duplicates,
  //           size <= base-domain-size?
  static void initialize(const std::vector<std::size_t>& indices, const std::string& name = "");

  // Use '_' (underscore) to be consistent with other domain classes.
  static std::size_t get_size() {
    assert(initialized_);
    return elements_.size();
  }
  static const std::string& get_name() {
    assert(initialized_);
    return name_;
  }
  static const std::vector<ElementType>& get_elements() {
    assert(initialized_);
    return elements_;
  }

private:
  static void commonInit(const std::string& name);

  static bool initialized_;
  static std::string name_;
  static std::vector<ElementType> elements_;
};

template <class BaseDomain>
bool ReducedDomain<BaseDomain>::initialized_ = false;
template <class BaseDomain>
std::string ReducedDomain<BaseDomain>::name_ = "";
template <class BaseDomain>
std::vector<typename ReducedDomain<BaseDomain>::ElementType> ReducedDomain<BaseDomain>::elements_;

template <class BaseDomain>
void ReducedDomain<BaseDomain>::initialize(const std::size_t n, const bool from_middle,
                                           const std::string& name) {
  commonInit(name);

  const std::size_t base_size = BaseDomain::dmn_size();

  if (from_middle) {
    if (n > (base_size + 1) / 2)
      throw std::out_of_range("Range to copy goes beyond end of base domain.");
  }
  else {
    if (n > base_size)
      throw std::out_of_range("Range to copy goes beyond end of base domain.");
  }

  elements_.assign(BaseDomain::get_elements().begin() + from_middle * base_size / 2,
                   BaseDomain::get_elements().begin() + from_middle * base_size / 2 + n);
}

template <class BaseDomain>
void ReducedDomain<BaseDomain>::initialize(const std::vector<std::size_t>& indices,
                                           const std::string& name) {
  commonInit(name);

  if (indices.size() > BaseDomain::dmn_size())
    throw std::out_of_range("More indices than elements in the base domain.");

  elements_.resize(indices.size());

  for (std::size_t i = 0; i < indices.size(); i++) {
    if (indices[i] >= BaseDomain::dmn_size())
      throw std::out_of_range("Indices must be < base-domain-size.");
    elements_[i] = BaseDomain::get_elements()[indices[i]];
  }
}

template <class BaseDomain>
void ReducedDomain<BaseDomain>::commonInit(const std::string& name) {
  if (initialized_)
    throw(std::logic_error("Domain has already been initialized."));

  initialized_ = true;

  // If no name ("") is passed, use BaseDomain_reduced.
  name_ = name == "" ? BaseDomain::get_name() + "_reduced" : name;
}

// Don't allow reduced domains of composed (variadic) domains.
template <class BaseDomain>
class ReducedDomain<dmn_variadic<BaseDomain>>;

}  // func
}  // dca

#endif  // COMP_LIBRARY_FUNCTION_LIBRARY_DOMAINS_SPECIAL_DOMAINS_REDUCED_DOMAIN_HPP
