// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides an alternative output mode which produces minimal amount of information about
// tests. Only failed assertions are disclosed.

#ifndef DCA_TESTING_MINIMALIST_PRINTER_HPP
#define DCA_TESTING_MINIMALIST_PRINTER_HPP

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

#endif  // DCA_TESTING_MINIMALIST_PRINTER_HPP
