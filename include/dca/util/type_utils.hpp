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
// This file provides utility functions that operate on types.
//
// TODO: Fix such that print_type.hpp can be removed.

#ifndef DCA_UTIL_TYPE_UTILS_HPP
#define DCA_UTIL_TYPE_UTILS_HPP

#include <complex>
#include <memory>
#include <string>
#include <type_traits>
#include <typeinfo>
#ifndef _MSC_VER
#include <cxxabi.h>
#endif

#ifdef DCA_HAVE_GPU
#include "dca/platform/dca_gpu_complex.h"
#endif

#include "dca/util/type_list.hpp"

namespace dca {
namespace func {
// dca::func::

// Forward declare these templates so we can use them in print functions.
// TODO: Move the domain print type functions to the domain classes.
template <typename Parameters>
class dmn_0;
template <typename... DomainList>
class dmn_variadic;
}  // namespace func

namespace util {
// dca::util::

// Throws an assertion if the two types do not match.
// Extends std::is_same<> by forcing the compiler to print the types in the error message which
// helps with debugging.
template <typename T1, typename T2>
struct assert_same {
  assert_same() {
    static_assert(std::is_same<T1, T2>::value, "Types must be equal.");
  }
  static_assert(std::is_same<T1, T2>::value, "Types must be equal.");
};

// Prints a type cleanly if possible.
// Reference:
// http://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c
template <class T>
std::string type_name() {
  using TR = typename std::remove_reference<T>::type;
  std::unique_ptr<char, void (*)(void*)> own(
#ifndef _MSC_VER
      abi::__cxa_demangle(typeid(TR).name(), nullptr, nullptr, nullptr),
#else
      nullptr,
#endif
      std::free);
  std::string r = own != nullptr ? own.get() : typeid(TR).name();
  if (std::is_const<TR>::value)
    r += " const";
  if (std::is_volatile<TR>::value)
    r += " volatile";
  if (std::is_lvalue_reference<T>::value)
    r += "&";
  else if (std::is_rvalue_reference<T>::value)
    r += "&&";
  return r;
}

// print_type
template <class T>
struct print_type_name {
  static void print() {
    print(std::cout);
  }

  static void print(std::ostream& ss) {
    ss << "\t" << T::get_name() << std::endl;
  }

  template <class stream_type>
  static void to_JSON(stream_type& ss) {
    ss << "\t" << T::get_name() << std::endl;
  }
};

// Basic print
// Displays the type of the template instantiation.
template <typename D>
struct print_type {
  static void print() {
    print(std::cout);
  }
  static void print(std::ostream& stream) {
    stream << "\t" << type_name<D>().c_str() << "\n";
  }
  static void to_JSON(std::ostream& stream) {
    stream << "\t" << type_name<D>().c_str();
  }
};

// dmn_0
// It is actually the same as basic, but provided for future customization.
template <typename Domain>
struct print_type<func::dmn_0<Domain>> {
  static void print() {
    print(std::cout);
  }
  static void print(std::ostream& stream) {
    stream << "\t" << type_name<func::dmn_0<Domain>>().c_str() << "\n";
  }
};

// dmn_variadic
// Prints out all subdomains recursively in left-to-right order via pack expansion.
template <typename... Domains>
struct print_type<func::dmn_variadic<Domains...>> {
  static void print() {
    print(std::cout);
  }

  // The std::initialize_list guarantees that the expressions that calculate its arguments
  // are executed in left-to-right order. It should have no overhead.
  // (http://florianjw.de/en/variadic_templates.html)
  // The void cast inside the initializer-list is a a precaution against potential
  // overloaded comma operators. Acutally it would be nicer to use lambda expressions
  // instead of the comma operator but a bug in gcc prevents us from doing this.
  static void print(std::ostream& s) {
    (void)std::initializer_list<int>{((void)print_type<Domains>::print(s), 0)...};
  }
};

// Typelist
// Prints out all subdomains recursively in left-to-right oder via pack expansion.
template <typename Domain, typename... Domains>
struct print_type<dca::util::Typelist<Domain, Domains...>> {
  static void print() {
    print(std::cout);
  }
  // Reuse dmn_variadic variant.
  static void print(std::ostream& s) {
    print_type<func::dmn_variadic<Domain, Domains...>>::print(s);
  }

  static void to_JSON(std::ostream& s) {
    s << "\"";
    print_type<Domain>::to_JSON(s);
    if (sizeof...(Domains) == 0) {
      s << "\"\n";
    }
    else {
      s << "\",\n";
      print_type<dca::util::Typelist<Domains...>>::to_JSON(s);
    }
  }
};

template <typename T>
struct IsComplex_t : public std::false_type {};
template <typename T>
struct IsComplex_t<std::complex<T>> : public std::true_type {};

#ifdef DCA_HAVE_GPU

template <>
struct IsComplex_t<cuComplex> : public std::true_type {};

template <>
struct IsComplex_t<cuDoubleComplex> : public std::true_type {};
#endif

template <typename T>
using IsComplex = std::enable_if_t<IsComplex_t<T>::value, bool>;
template <typename T>
using IsReal = std::enable_if_t<std::is_floating_point<T>::value, bool>;

template <typename T, typename = bool>
struct RealAlias_impl {};

template <typename T>
struct RealAlias_impl<T, IsReal<T>> {
  using value_type = T;
};

#ifdef DCA_HAVE_GPU
template <>
struct RealAlias_impl<cuComplex, bool> {
  using value_type = float;
};

template <>
struct RealAlias_impl<cuDoubleComplex, bool> {
  using value_type = double;
};
#endif

template <typename T>
struct RealAlias_impl<T, IsComplex<T>> {
  using value_type = typename T::value_type;
};

/** If you have a function templated on a value that can be real or complex
 *   and you need to get the base Real type if its complex or just the real.
 *
 *  If you try to do this on anything but a fp or a std::complex<fp> you will
 *  get a compilation error.
 */
template <typename T>
using RealAlias = typename RealAlias_impl<T>::value_type;

template <typename T, typename = bool>
struct ComplexAlias_impl {};

template <typename T>
struct ComplexAlias_impl<T, IsReal<T>> {
  using value_type = std::complex<T>;
};

template <typename T>
struct ComplexAlias_impl<T, IsComplex<T>> {
  using value_type = T;
};

template <typename T>
using ComplexAlias = typename ComplexAlias_impl<T>::value_type;

template <typename REAL, bool complex>
struct ScalarSelect {
  using type = REAL;
};

template <typename REAL>
struct ScalarSelect<REAL, false> {
  using type = REAL;
};

template <typename REAL>
struct ScalarSelect<REAL, true> {
  using type = std::complex<REAL>;
};

template <typename T, typename = bool>
struct TheOne {
  static constexpr T value = 1;
};

template <typename T>
struct TheOne<T, IsComplex<T>> {
  static constexpr T value = {1.0, 0.0};
};

template <typename T, typename = bool>
struct TheZero {
  static constexpr T value = 0;
};

template <typename T>
struct TheZero<T, IsComplex<T>> {
  static constexpr T value = {0.0, 0.0};
};

template <typename ARRAY, std::size_t... SIZE>
auto Array2Tuple_impl(const ARRAY& a, std::index_sequence<SIZE...>) {
  return std::make_tuple(a[SIZE]...);
}

template <typename T, std::size_t N, typename Indices = std::make_index_sequence<N>>
auto Array2Tuple(const std::array<T, N>& a) {
  return Array2Tuple_impl(a, Indices{});
}

}  // namespace util
}  // namespace dca

#endif  // DCA_UTIL_TYPE_UTILS_HPP
