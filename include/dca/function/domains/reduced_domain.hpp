// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements a class to store a subset of a previously defined one-dimensional domain.
// This subset can either be a set of single elements or a range of the original domain.
// The domain is created by a type definition followed by a call to an initialize function.
//
// INTERNAL: This class is probably good enough for its purpose, but it doesn't seem to be
//           'complete'.
//
// TODO: This class can only be used with dmn_variadic when it's wrapped with dmn_0.

#ifndef DCA_FUNCTION_DOMAINS_REDUCED_DOMAIN_HPP
#define DCA_FUNCTION_DOMAINS_REDUCED_DOMAIN_HPP

#include <algorithm>  // for std::is_sorted
#include <cassert>
#include <cstdlib>  // for std::size_t
#include <stdexcept>
#include <string>
#include <vector>

#include "dca/function/domains/dmn_variadic.hpp"

namespace dca {
namespace func {
// dca::func::

template <typename BaseDomain>
class ReducedDomain {
public:
  using ElementType = typename BaseDomain::element_type;
  using element_type = ElementType;  // For putting ReducedDomain in dmn_0.

  // Creates a domain that contains the elements in the range [first, last) of the base domain.
  static void initialize(const std::size_t first, const std::size_t last,
                         const std::string& name = "");

  // Creates a domain that only contains the elements of the base domain with the given indices.
  // Preconditions: - 0 <= indices[i] < base-domain-size, for 0 <= i < indices.size().
  //                - indices is sorted and doesn't contain duplicates.
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
void ReducedDomain<BaseDomain>::initialize(const std::size_t first, const std::size_t last,
                                           const std::string& name) {
  commonInit(name);

  const std::size_t base_size = BaseDomain::dmn_size();

  if (first >= base_size)
    throw std::out_of_range("First index must be < base-domain-size.");
  if (last > base_size)
    throw std::out_of_range("Last index must be <= base-domain-size.");
  if (first > last)
    throw std::logic_error("Last index must be >= first index.");

  elements_.assign(BaseDomain::get_elements().begin() + first,
                   BaseDomain::get_elements().begin() + last);
}

template <class BaseDomain>
void ReducedDomain<BaseDomain>::initialize(const std::vector<std::size_t>& indices,
                                           const std::string& name) {
  commonInit(name);

  if (!std::is_sorted(indices.begin(), indices.end()))
    throw std::logic_error("Indices must be sorted.");

  elements_.resize(indices.size());

  for (std::size_t i = 0; i < indices.size(); i++) {
    if (indices[i] >= BaseDomain::dmn_size())
      throw std::out_of_range("Indices must be < base-domain-size.");
    if (i > 0 && indices[i - 1] == indices[i])
      throw std::logic_error("Vector with indices contains duplicates.");

    elements_[i] = BaseDomain::get_elements()[indices[i]];
  }
}

template <class BaseDomain>
void ReducedDomain<BaseDomain>::commonInit(const std::string& name) {
  if (initialized_)
    throw std::logic_error("Domain has already been initialized.");

  initialized_ = true;

  // If no name ("") is passed, use BaseDomain_reduced.
  name_ = name == "" ? BaseDomain::get_name() + "_reduced" : name;
}

// Don't allow reduced domains of composed (variadic) domains.
template <class BaseDomain>
class ReducedDomain<dmn_variadic<BaseDomain>>;

}  // func
}  // dca

#endif  // DCA_FUNCTION_DOMAINS_REDUCED_DOMAIN_HPP
