// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides an alternative output mode which produces minimal amount of information about
// tests. Only failed assertions are disclosed.

#ifndef DCA_TEST_COMMON_MINIMALIST_PRINTER_HPP
#define DCA_TEST_COMMON_MINIMALIST_PRINTER_HPP

#include "gtest/gtest.h"

namespace dca {
namespace testing {
// dca::testing::

class MinimalistPrinter : public ::testing::EmptyTestEventListener {
private:
  // Called after a failed assertion or a SUCCEED() invocation.
  virtual void OnTestPartResult(const ::testing::TestPartResult& test_part_result) {
    if (test_part_result.failed()) {
      fprintf(stdout, "%s in %s:%d\n%s\n", "*** Failure", test_part_result.file_name(),
              test_part_result.line_number(), test_part_result.summary());
      fflush(stdout);
    }
  }
};

}  // testing
}  // dca

#endif  // DCA_TEST_COMMON_MINIMALIST_PRINTER_HPP
