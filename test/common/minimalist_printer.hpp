// Author: Urs Haehner (haehneru@itp.phys.ethz.ch), ETH Zurich

#ifndef MINIMALIST_PRINTER_HPP
#define MINIMALIST_PRINTER_HPP

// MinimalistPrinter class
// Provides alternative output mode which produces minimal amount of information
// about tests. Only failed assertions are disclosed.
class MinimalistPrinter : public ::testing::EmptyTestEventListener {
private:
  // Called after a failed assertion or a SUCCEED() invocation.
  virtual void OnTestPartResult(const ::testing::TestPartResult& test_part_result) {
    if (test_part_result.failed()) {
      fprintf(stdout,
              "%s in %s:%d\n%s\n",
              "*** Failure",
              test_part_result.file_name(),
              test_part_result.line_number(),
              test_part_result.summary());
      fflush(stdout);
    }
  }
};

#endif  // MINIMALIST_PRINTER_HPP
