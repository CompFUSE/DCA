// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This class provides an alternative output mode which produces minimal amount of information about
// tests. Only failed assertions are disclosed.

#ifndef DCA_TESTING_TYPE_TESTING_HPP
#define DCA_TESTING_TYPE_TESTING_HPP

#include <type_traits>

namespace dca {
namespace testing {

template<typename...>
using void_t = void;

template<typename, typename, typename = void>
struct can_compare : std::false_type {};

template<typename A, typename B>
struct can_compare<A, B,
    void_t<decltype(std::declval<A>() == std::declval<B>())>
                   > : std::true_type {};

}
}
#endif
