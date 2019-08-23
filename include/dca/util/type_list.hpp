// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: John Biddiscombe (john.biddiscombe@cscs.ch)
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

// mp_rename: A<...> -> B<...>
template <class A, template <class...> class B>
struct mp_rename_impl;

template <template <class...> class A, class... T, template <class...> class B>
struct mp_rename_impl<A<T...>, B> {
  using type = B<T...>;
};

template <class A, template <class...> class B>
using mp_rename = typename mp_rename_impl<A, B>::type;

// mp_size
template <class L>
struct mp_size_impl;

template <template <class...> class L, class... T>
struct mp_size_impl<L<T...>> {
  using type = std::integral_constant<std::size_t, sizeof...(T)>;
};

template <class L>
using mp_size = typename mp_size_impl<L>::type;

// mp_plus
template <class... T>
struct mp_plus_impl;

template <class... T>
using mp_plus = typename mp_plus_impl<T...>::type;

template <>
struct mp_plus_impl<> {
  using type = std::integral_constant<int, 0>;
};

template <class T1, class... T>
struct mp_plus_impl<T1, T...> {
  static constexpr auto _v = T1::value + mp_plus<T...>::value;
  using type = std::integral_constant<typename std::remove_const<decltype(_v)>::type, _v>;
};

// mp_count
template <class L, class V>
struct mp_count_impl;

template <template <class...> class L, class... T, class V>
struct mp_count_impl<L<T...>, V> {
  using type = mp_plus<std::is_same<T, V>...>;
};

template <class L, class V>
using mp_count = typename mp_count_impl<L, V>::type;

// mp_append
template <typename T, typename... TL>
struct mp_append;

template <typename... Ts>
struct mp_append<mp_list<Ts...>> {
  typedef mp_list<Ts...> type;
};

template <typename T, typename... Ts>
struct mp_append<mp_list<Ts...>, T> {
  typedef mp_list<Ts..., T> type;
};

template <typename... Ts1, typename... Ts2>
struct mp_append<mp_list<Ts1...>, mp_list<Ts2...>> {
  typedef mp_list<Ts1..., Ts2...> type;
};

template <typename... Ts1, typename... Ts2, typename... Ts>
struct mp_append<mp_list<Ts1...>, mp_list<Ts2...>, Ts...> {
  typedef typename mp_append<mp_list<Ts1..., Ts2...>, Ts...>::type type;
};

// mp_prepend
template <typename T, typename... TL>
struct mp_prepend;

template <typename T, typename... Ts>
struct mp_prepend<mp_list<Ts...>, T> {
  typedef mp_list<T, Ts...> type;
};

template <typename... Ts1, typename... Ts2>
struct mp_prepend<mp_list<Ts1...>, mp_list<Ts2...>> {
  typedef mp_list<Ts2..., Ts1...> type;
};

// mp_element
// Get the n'th type from a typelist/tuple efficiently without recursion.
// Reference: True Story: Efficient Packing
//            http://talesofcpp.fusionfenix.com/post-22/true-story-efficient-packing
template <std::size_t I, typename T>
struct _indexed {
  using type = T;
};

template <typename Is, typename... Ts>
struct _indexer;

template <std::size_t... Is, typename... Ts>
struct _indexer<std::index_sequence<Is...>, Ts...> : _indexed<Is, Ts>... {};

template <std::size_t I, typename... Ts>
struct _at_index {
  template <typename T>
  static _indexed<I, T> _select(_indexed<I, T>);

  using _impl = _indexer<std::index_sequence_for<Ts...>, Ts...>;
  using type = typename decltype(_select(_impl{}))::type;
};

// Was tuple_element, but we want typelist_element.
template <std::size_t I, typename Tuple>
struct mp_element;

template <std::size_t I, typename... Ts>
struct mp_element<I, mp_list<Ts...>> : _at_index<I, Ts...> {};

// mp_index_of
// Search a typelist for a first occurrence of the type T.
// Implementation: has index as a template parameter
template <size_t idx, typename T, class List>
struct mp_index_of_impl;

template <size_t idx, typename T>  // The type T is not in the list.
struct mp_index_of_impl<idx, T, mp_list<>> {
  using type = std::integral_constant<int, -1>;
};

template <size_t idx, typename T, typename... Ts>  // The type is found.
struct mp_index_of_impl<idx, T, mp_list<T, Ts...>> {
  using type = std::integral_constant<int, idx>;
};

template <size_t idx, typename T, typename H, typename... Ts>  // Recursion.
struct mp_index_of_impl<idx, T, mp_list<H, Ts...>> {
  using type = typename mp_index_of_impl<idx + 1, T, mp_list<Ts...>>::type;
};

// Wrapping to supply initial index 0.
template <typename T, class List>
struct mp_index_of {
  static constexpr int value = -1;
};

// Specializing for idx >= 0.
template <typename T, typename... Ts>
struct mp_index_of<T, mp_list<Ts...>> {
  using type = typename mp_index_of_impl<0, T, mp_list<Ts...>>::type;
  using value_type = typename type::value_type;
  static constexpr value_type value = type::value;
};

// Specializing for idx >= 0.
template <typename T, typename... Ts>
struct mp_index_of<mp_list<Ts...>, T> {
  // static_assert(false, "Parameter ordering incorrect");
};

// mp_swap
// Swap first element out of a typelist.
template <typename TList, typename T1, typename T2>
struct mp_swap {};

template <typename T1, typename T2>
struct mp_swap<mp_list<>, T1, T2> {
  typedef mp_list<> type;
};

template <typename T1, typename T2, typename... Ts>
struct mp_swap<mp_list<T1, Ts...>, T1, T2> {
  typedef mp_list<T2, Ts...> type;
};

template <typename T0, typename... Ts, typename T1, typename T2>
struct mp_swap<mp_list<T0, Ts...>, T1, T2> {
  typedef typename mp_prepend<typename mp_swap<mp_list<Ts...>, T1, T2>::type, T0>::type type;
};

template <typename T1, typename Ts>
constexpr bool contained() {
  return mp_index_of<T1, Ts>::value != -1;
}

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

}  // namespace util
}  // namespace dca

#endif  // DCA_UTIL_TYPE_LIST_HPP
