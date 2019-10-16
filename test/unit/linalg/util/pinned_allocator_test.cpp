// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file tests the pinned allocator

#include "dca/linalg/util/allocators/pinned_allocator.hpp"
#include "dca/testing/type_testing.hpp"
#include <stdexcept>
#include <cuda_runtime.h>
#include "gtest/gtest.h"

TEST(PinnedAllocatorTest, Allocator) {
  using dca::linalg::util::PinnedAllocator;
  std::vector<int, PinnedAllocator<int>> pinned_array1(20,1.0);
  std::vector<int, PinnedAllocator<int>> pinned_array2;
  pinned_array2 = std::move(pinned_array1);

  PinnedAllocator<double> dalloc;
  PinnedAllocator<double> d2alloc;
  PinnedAllocator<int> ialloc;
  double* mymem = dalloc.allocate(64);
  std::allocator<double> sdalloc;
  EXPECT_TRUE( dalloc == d2alloc );
  EXPECT_TRUE( dalloc == ialloc );

  if(dalloc == d2alloc)
    dalloc.deallocate(mymem);

  // Test that compile time logic works as expected
  static_assert(! dca::testing::can_compare<decltype(dalloc), decltype(sdalloc)>::value, "PinnedAllocator and std::allocator should not be comparable");
  static_assert( dca::testing::can_compare<decltype(dalloc), decltype(d2alloc)>::value, "PinnedAllocator's should be comparable");
  static_assert( dca::testing::can_compare<decltype(dalloc), decltype(ialloc)>::value, "PinnedAllocator's should be comparable");
  
}
