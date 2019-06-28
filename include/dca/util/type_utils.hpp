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
}  // func

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

// Determine if a type is complex or not.
template <class T>
struct IsComplex {
  constexpr static bool value = 0;
};
template <class T>
struct IsComplex<std::complex<T>> {
  constexpr static bool value = 1;
};

}  // util
}  // dca

#endif  // DCA_UTIL_TYPE_UTILS_HPP
