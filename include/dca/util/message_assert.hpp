// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W.  Doak (doakpw@ornl.gov)
//

#ifndef DCA_UTIL_MESSAGE_ASSERT_HPP
#define DCA_UTIL_MESSAGE_ASSERT_HPP

#include <cassert>
#include <iostream>

#define MESSAGE_ASSERT(condition, message)                                                        \
  do {                                                                                            \
    if (!(condition)) {                                                                           \
      std::cerr << "Assertion failed: " << #condition << "\nMessage: " << (message) << std::endl; \
      assert(false);                                                                              \
    }                                                                                             \
  } while (0)

#endif
