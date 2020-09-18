// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: John Biddiscombe (john.biddiscombe@cscs.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a type list and type list operations.
//
// References: - Simple C++11 metaprogramming, Peter Dimov
//               http://pdimov.com/cpp2/simple_cxx11_metaprogramming.html
//               Interesting discussion of type transformations and examples.
//               mp_xxxx metafunctions are taken from this source.
//             - TypeLists and a TypeList Toolbox via Variadic Templates
//               http://www.codeproject.com/Articles/1077852/TypeLists-and-a-TypeList-Toolbox-via-Variadic-Temp

#ifndef DCA_UTIL_TYPE_LIST_HPP
#define DCA_UTIL_TYPE_LIST_HPP

#include <iostream>
#include <type_traits>
#include <utility>

namespace dca {
namespace util {
// dca::util::

// mp_list: generic type list
template <class... T>
struct mp_list {};

// mp_size
template <class L>
struct mp_size;

template <template <class...> class L, class... T>
struct mp_size<L<T...>> {
  constexpr static int value = sizeof...(T);
};

// mp_count
template <class L, class V>
struct mp_count;

template <template <class...> class L, class... T, class V>
struct mp_count<L<T...>, V> {
  constexpr static int value = (0 + ... + std::is_same_v<T, V>);
};

// mp_append
template <typename T, typename... TL>
struct mp_append;

template <typename... Ts>
struct mp_append<mp_list<Ts...>> {
  typedef mp_list<Ts...> type;
};

template <typename T, typename... Ts>
struct mp_append<mp_list<Ts...>, T> {
  using type = mp_list<Ts..., T>;
};

template <typename... Ts1, typename... Ts2>
struct mp_append<mp_list<Ts1...>, mp_list<Ts2...>> {
  using type = mp_list<Ts1..., Ts2...>;
};

template <typename... Ts1, typename... Ts2, typename... Ts>
struct mp_append<mp_list<Ts1...>, mp_list<Ts2...>, Ts...> {
  using type = typename mp_append<mp_list<Ts1..., Ts2...>, Ts...>::type;
};

// mp_prepend
template <typename T, typename... TL>
struct mp_prepend;

template <typename T, typename... Ts>
struct mp_prepend<mp_list<Ts...>, T> {
  using type = mp_list<T, Ts...>;
};

template <typename... Ts1, typename... Ts2>
struct mp_prepend<mp_list<Ts1...>, mp_list<Ts2...>> {
  using type = mp_list<Ts2..., Ts1...>;
};

// mp_element
// Was tuple_element, but we want typelist_element.
template <std::size_t I, typename Tuple>
struct mp_element;

template <std::size_t I, typename... Ts>
struct mp_element<I, mp_list<Ts...>> {
  using type = std::tuple_element_t<I, std::tuple<Ts...>>;
};

// mp_index_of
// Search a typelist for a first occurrence of the type T. Value is -1 if the type is not found.
template <typename V, typename... Ts>
struct mp_index_of;

template <typename V, typename... Ts>
struct mp_index_of<V, mp_list<V, Ts...>> {  // Index found.
  static constexpr int value = 0;
};

template <typename V>
struct mp_index_of<V, mp_list<>> {  // Index not found.
  static constexpr int value = -1;
};

template <typename V, typename T1, typename... Ts>
struct mp_index_of<V, mp_list<T1, Ts...>> {  // Recursion.
  static constexpr int next = mp_index_of<V, mp_list<Ts...>>::value;
  static constexpr int value = next == -1 ? -1 : next + 1;
};

// mp_swap
// Swap first element out of a typelist.
template <typename TList, typename T1, typename T2>
struct mp_swap {};

template <typename T1, typename T2>
struct mp_swap<mp_list<>, T1, T2> {
  using type = mp_list<>;
};

template <typename T1, typename T2, typename... Ts>
struct mp_swap<mp_list<T1, Ts...>, T1, T2> {
  using type = mp_list<T2, Ts...>;
};

template <typename T0, typename... Ts, typename T1, typename T2>
struct mp_swap<mp_list<T0, Ts...>, T1, T2> {
  using type = typename mp_prepend<typename mp_swap<mp_list<Ts...>, T1, T2>::type, T0>::type;
};

template <typename T1, typename Domain>
constexpr bool contained() {
  return mp_index_of<T1, typename Domain::this_type>::value != -1;
}

// mp_sublist
// Creates a sublist with the first n elements of a given list.
template <unsigned n, typename T1, typename... Ts>
struct mp_sublist {
  using next = typename mp_sublist<n - 1, Ts...>::type;

  using type = typename mp_prepend<next, T1>::type;
};
template <typename T1, typename... Ts>
struct mp_sublist<0, T1, Ts...> {
  using type = mp_list<>;
};

template <unsigned id, typename... Ts>
struct mp_sublist<id, mp_list<Ts...>> {
  using type = typename mp_sublist<id, Ts...>::type;
};

// Create aliases for backward naming compatibility with old typelist and typelist operations.
template <class... Ts>
using Typelist = mp_list<Ts...>;

template <typename L>
using Length = mp_size<L>;

template <class L, class V>
using NumberOf = mp_count<L, V>;

template <int I, typename Ts>
using TypeAt = mp_element<I, Ts>;

template <typename T, typename Ts>
using IndexOf = mp_index_of<T, Ts>;

template <typename T1, typename T2>
using Append = mp_append<T1, T2>;

template <typename T1, typename T2>
using Prepend = mp_prepend<T1, T2>;

template <typename T1, typename T2, typename T3>
using Swap = mp_swap<T1, T2, T3>;

template <unsigned n, typename... Ts>
using Sublist = typename mp_sublist<n, Ts...>::type;

}  // namespace util
}  // namespace dca

#endif  // DCA_UTIL_TYPE_LIST_HPP
