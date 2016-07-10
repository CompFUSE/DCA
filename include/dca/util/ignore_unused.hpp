// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides a function to silence an unused parameter/variable compiler warning.
// Reference: Stack Overflow http://stackoverflow.com/questions/15763937/unused-parameter-in-c11

#ifndef DCA_UTIL_IGNORE_UNUSED_HPP
#define DCA_UTIL_IGNORE_UNUSED_HPP

namespace dca {
namespace util {
// dca::util::

template <class... T>
void ignoreUnused(T&&...) {}

}  // util
}  // dca

#endif  // DCA_UTIL_IGNORE_UNUSED_HPP
