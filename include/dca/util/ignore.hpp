// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         John Biddiscombe (john.biddiscombe@cscs.ch)
//
// This file provides a utility function that takes a variadic list of arguments and returns
// nothing.
// We provide this function with different names for code readability.

#ifndef DCA_UTIL_IGNORE_HPP
#define DCA_UTIL_IGNORE_HPP

namespace dca {
namespace util {
// dca::util::

// Silences an unused parameter/variable compiler warning.
// Reference: Stack Overflow http://stackoverflow.com/questions/15763937/unused-parameter-in-c11
template <typename... Ts>
void ignoreUnused(Ts&&...) {}

// Allows to call arbitrary functions with a variadic pack expansion.
template <typename... Ts>
void ignoreReturnValues(Ts&&...) {}

}  // util
}  // dca

#endif  // DCA_UTIL_IGNORE_HPP
