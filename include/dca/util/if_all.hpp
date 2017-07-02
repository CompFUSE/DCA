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
// if_all<b1, b2, ...>::value is true only if all template arguments are true, otherwise false.

#ifndef DCA_UTIL_IF_ALL_HPP
#define DCA_UTIL_IF_ALL_HPP

namespace dca {
namespace util {
// dca::util::

template <bool b1, bool... bs>
struct if_all {
  constexpr static bool value = b1 && if_all<bs...>::value;
};

template <bool b>
struct if_all<b> {
  constexpr static bool value = b;
};

}  // util
}  // dca

#endif  // DCA_UTIL_IF_ALL_HPP
