// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// No-change test for CT-INT posix wrapper.

using Scalar = double;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca


#include "stdthread_ctint_cluster_solver_test.hpp"

TEST(PosixCtintClusterSolverTest, Shared) {
  performBaselineTest("stdthread_ctint_test_shared_input.json",
                      "stdthread_ctint_test_shared_baseline.hdf5");
}
