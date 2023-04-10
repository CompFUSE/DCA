// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific
// publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)
//

#ifndef DCA_TEST_MOCK_MCCONFIG_HPP
#define DCA_TEST_MOCK_MCCONFIG_HPP

#include "dca/util/type_utils.hpp"
#include "dca/config/haves_defines.hpp"
#ifdef DCA_HAVE_GPU
#include "dca/linalg/util/allocators/device_allocator.hpp"
#include "dca/linalg/util/allocators/managed_allocator.hpp"
#endif  // DCA_HAVE_GPU


namespace dca {
namespace config {

// Mocking global build config for tests
template <typename Scalar>
struct MockMcOptions {
  using MC_REAL = dca::util::RealAlias<Scalar>;
  using TPAccumulationPrecision = dca::util::RealAlias<Scalar>;
  static constexpr bool memory_savings = false;
#ifdef DCA_HAVE_GPU
  template <typename T>
  using TpAllocator = dca::linalg::util::DeviceAllocator<T>;
#endif  // DCA_HAVE_GPU
};

}  // namespace config
}  // namespace dca

#endif
