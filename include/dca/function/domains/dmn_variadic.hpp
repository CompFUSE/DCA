// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: John Biddiscombe (john.biddiscombe@cscs.ch)
//
// Variadic domain class with a list of sub-domains.

#ifndef DCA_FUNCTION_DOMAINS_DMN_VARIADIC_HPP
#define DCA_FUNCTION_DOMAINS_DMN_VARIADIC_HPP

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <string>
#include <tuple>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "dca/function/domains/domain.hpp"
#include "dca/util/ignore.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace func {
// dca::func::

template <typename... domain_list>
class dmn_variadic : public domain {
public:
  typedef typename dca::util::mp_append<typename domain_list::this_type...>::type this_type;

  template <int Index>
  using domain_typelist =
      typename std::tuple_element<Index, std::tuple<typename domain_list::this_type...>>::type;

  // Constructor, responsible for initializing all domain indexing arrays.
  dmn_variadic();

  // Returns the size allocated by this domain.
  // Includes all space for all subdomains.
  static int& dmn_size() {
    static int size = -1;
    // std::cout << "Returning domain size " << size << std::endl << std::endl;
    return size;
  }

  // Resets the domain back to the original state at creation time.
  void reset();

  // Indexing operator: accesses elements of the domain.
  // If sizeof(args) < sizeof(domains): error.
  // If sizeof(args) == sizeof(domains): calls index via branch domains.
  // If sizeof(args) == sizeof(leaf domains): calls index via leaf domains.
  template <typename... Args>
  std::size_t operator()(Args&&... args) const;

protected:
  // Indexing operator: accesses elements of the domain via branches.
  // index_lookup is overloaded on std::integral_constant<bool, true>::type so that if
  // sizeof...(Args) == sizeof...(domain_list) then this is called.
  template <typename... Args>
  std::size_t index_lookup(std::integral_constant<bool, true>::type, int branch_i0,
                           Args... args) const;

  // Indexing operator: accesses elements of the domain via leaf domains.
  // index_lookup is overloaded on std::integral_constant<bool, true>::type so that if
  // sizeof...(Args) == sizeof...(domain_list) then this is called.
  template <typename... Args>
  std::size_t index_lookup(std::integral_constant<bool, false>::type, int leaf_i0, Args... args) const;

  // Gets the branch domain sizes for each of the domain template arguments.
  // Returns a vector of the sizes.
  template <typename Tuple, std::size_t... DomainIndex>
  std::vector<std::size_t> init_branch_domain_sizes(Tuple& t, std::index_sequence<DomainIndex...>);

  // Gets the number of leaf domains for each of the domain template arguments.
  // Returns a vector of the counts.
  template <typename Tuple, std::size_t... DomainIndex>
  std::vector<std::size_t> init_leaf_domain_sizes(Tuple& domaintuple,
                                                  std::index_sequence<DomainIndex...>);

protected:
  std::tuple<domain_list...> domains;
};

template <typename... domain_list>
template <typename Tuple, std::size_t... DomainIndex>
std::vector<std::size_t> dmn_variadic<domain_list...>::init_branch_domain_sizes(
    Tuple& domaintuple, std::index_sequence<DomainIndex...>) {
  return {(std::get<DomainIndex>(domaintuple).get_size())...};
}

template <typename... domain_list>
template <typename Tuple, std::size_t... DomainIndex>
std::vector<std::size_t> dmn_variadic<domain_list...>::init_leaf_domain_sizes(
    Tuple& domaintuple, std::index_sequence<DomainIndex...>) {
  return {(std::get<DomainIndex>(domaintuple).get_Nb_leaf_domains())...};
}

namespace detail {

// Invokes a function for each element of a tuple.
template <typename T, typename F, std::size_t... Indices>
void for_each(T&& t, F&& f, std::index_sequence<Indices...>) {
  // We need to create this unused object to prevent a runtime error with GCC.
  auto l = {(f(std::get<Indices>(t)), 0)...};
  dca::util::ignoreUnused(l);
}

template <typename... Ts, typename F>
void for_each_in_tuple(std::tuple<Ts...>& t, F&& f) {
  for_each(t, std::forward<F>(f), std::make_index_sequence<sizeof...(Ts)>{});
}

// Modifies an index sequence, offset/length.
// For indexing operations we multiply arguments by branch/leaf step sizes.
template <std::size_t O, std::size_t... Is>
std::index_sequence<(O + Is)...> add_offset(std::index_sequence<Is...>) {
  return {};
}

template <std::size_t O, std::size_t N>
auto make_index_sequence_with_offset() {
  return add_offset<O>(std::make_index_sequence<N>{});
}

// Index multiplication
template <typename T>
T sum(T v) {
  return v;
}

template <typename T, typename... Args>
T sum(T first, Args... args) {
  return first + sum(args...);
}

template <typename... Args, std::size_t... Is>
std::size_t multiply_offsets(const std::vector<std::size_t>& multipliers,
                             std::index_sequence<Is...>, Args&&... offsets) {
  return sum((offsets * (multipliers[Is]))...);
}

}  // namespace detail

template <typename... domain_list>
dmn_variadic<domain_list...>::dmn_variadic() : domain() {
  dmn_variadic<domain_list...>::reset();
}

template <typename... domain_list>
void dmn_variadic<domain_list...>::reset() {
  domain::reset();

  detail::for_each_in_tuple(domains, [](auto& d) { d.reset(); });

  // Create an index sequence that indexes the domains we are templated on.
  std::index_sequence_for<domain_list...> indices;

  // Initialize a vector from the size of each top level domain (branch domain).
  branch_domain_sizes = init_branch_domain_sizes(domains, indices);

  // std::cout << "Creating " << __PRETTY_FUNCTION__ << " " << branch_domain_sizes.size() <<
  // std::endl
  //           << "domain sizes : ";
  // std::copy(branch_domain_sizes.begin(), branch_domain_sizes.end(),
  //           std::ostream_iterator<std::size_t>(std::cout, ","));
  // std::cout << std::endl;

  branch_domain_steps.resize(branch_domain_sizes.size(), 1);
  for (size_t i = 0; i < branch_domain_sizes.size(); i++) {
    for (size_t j = 0; j < i; j++) {
      branch_domain_steps[i] *= branch_domain_sizes[j];
    }
  }

  // std::cout << "Steps ";
  // std::copy(branch_domain_steps.begin(), branch_domain_steps.end(),
  //           std::ostream_iterator<std::size_t>(std::cout, ","));
  // std::cout << std::endl;

  // Generate the leaf domain sizes from each sub domain.
  auto leaf_stuff = init_leaf_domain_sizes(domains, indices);
  // std::cout << "Leaf stuff ";
  // std::copy(leaf_stuff.begin(), leaf_stuff.end(), std::ostream_iterator<std::size_t>(std::cout,
  // ",")); std::cout << std::endl;

  detail::for_each_in_tuple(domains, [this](auto& d) {
    std::vector<std::size_t> temp = d.get_leaf_domain_sizes();
    std::copy(std::begin(temp), std::end(temp), std::back_inserter(leaf_domain_sizes));
  });

  if (leaf_domain_sizes.back() == 0) {
    // throw std::runtime_error("Domain size zero");
  }

  leaf_domain_steps.resize(leaf_domain_sizes.size(), 1);
  for (size_t i = 0; i < leaf_domain_sizes.size(); i++)
    for (size_t j = 0; j < i; j++)
      leaf_domain_steps[i] *= leaf_domain_sizes[j];

  // std::cout << "leaf step size ";
  // std::copy(leaf_domain_steps.begin(), leaf_domain_steps.end(),
  //           std::ostream_iterator<std::size_t>(std::cout, ","));
  // std::cout << std::endl;

  size = 1;

  for (int i = 0; i < sizeof...(domain_list); i++)
    size *= branch_domain_sizes[i];

  dmn_size() = size;
  // std::cout << "Domain size is " << size << std::endl;
}

template <typename... domain_list>
template <typename... Args>
std::size_t dmn_variadic<domain_list...>::operator()(Args&&... args) const {
  static_assert(sizeof...(Args) >= sizeof...(domain_list), "not enough args");
  return index_lookup(std::integral_constant<bool, (sizeof...(Args) == sizeof...(domain_list))>(),
                      std::forward<Args>(args)...);
}

template <typename... Args, std::size_t... Is>
void check_indices(const char* /*msg*/, const std::vector<std::size_t>& sizes,
                   std::index_sequence<Is...>, Args&&... indices) {
  if (std::min({(sizes[Is] - indices)...}) <= 0) {
    dca::util::ignoreReturnValues(
        (std::cerr << "size " << sizes[Is] << " index " << indices << " ")...);
    std::cerr << " : Index too big error" << std::endl;
    std::copy(sizes.begin(), sizes.end(), std::ostream_iterator<std::size_t>(std::cerr, ","));
    throw std::runtime_error("Index too big error");
  }
  if (std::min({indices...}) < 0) {
    std::cerr << "Index too small error" << std::endl;
    throw std::runtime_error("Index too small error");
  }
}

template <typename... domain_list>
template <typename... Args>
std::size_t dmn_variadic<domain_list...>::index_lookup(std::integral_constant<bool, true>::type,
                                                       int branch_i0, Args... branch_indices) const {
  static_assert(sizeof...(Args) + 1 == sizeof...(domain_list), "not enough args");

  // Create an index sequence starting from 1, with length sizeof...(args)-1.
  auto seq = detail::make_index_sequence_with_offset<1, sizeof...(Args)>();

#ifndef NDEBUG
  auto seq2 = std::make_index_sequence<sizeof...(Args) + 1>{};
  check_indices("branch ", branch_domain_sizes, seq2, branch_i0,
                std::forward<Args>(branch_indices)...);
#endif  // NDEBUG

  std::size_t N = branch_i0 + detail::multiply_offsets(branch_domain_steps, seq,
                                                       std::forward<Args>(branch_indices)...);
  return N;
}

template <typename... domain_list>
template <typename... Args>
std::size_t dmn_variadic<domain_list...>::index_lookup(std::integral_constant<bool, false>::type,
                                                       int leaf_i0, Args... leaf_indices) const {
  // Create an index sequence starting from 1, with length sizeof...(args)-1.
  auto seq = detail::make_index_sequence_with_offset<1, sizeof...(Args)>();

#ifndef NDEBUG
  auto seq2 = std::make_index_sequence<sizeof...(Args) + 1>{};
  check_indices("leaf", leaf_domain_sizes, seq2, leaf_i0, std::forward<Args>(leaf_indices)...);
#endif  // NDEBUG

  std::size_t N = leaf_i0 + detail::multiply_offsets(leaf_domain_steps, seq,
                                                     std::forward<Args>(leaf_indices)...);
  return N;
}

}  // namespace func
}  // namespace dca

#endif  // DCA_FUNCTION_DOMAINS_DMN_VARIADIC_HPP
