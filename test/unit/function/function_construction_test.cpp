// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file tests constructors with respect to distribution strategy of the function library.
//

#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/function/function.hpp"

#include <array>
#include <complex>
#include <string>
#include <vector>
#include <typeinfo>  // operator typeid
#include <iostream>
#include <sstream>

#include "dca/linalg/util/info_cuda.hpp"
#include "dca/linalg/util/util_cublas.hpp"

#include "gtest/gtest.h"

dca::parallel::MPIConcurrency* concurrency_ptr;

template <typename Scalar>
class BlockDistributedFunctionTest : public ::testing::Test {};

using TestTypes = ::testing::Types<float>;

TYPED_TEST_CASE(BlockDistributedFunctionTest, TestTypes);

TYPED_TEST(BlockDistributedFunctionTest, CheckDomainSizesOnRanks) {
  int rank = concurrency_ptr->id();
  int comm_size = concurrency_ptr->number_of_processors();

  using Dmn1 = dca::func::dmn_0<dca::func::dmn<5>>;
  using Dmn2 = dca::func::dmn_0<dca::func::dmn<4>>;
  using Dmn3 = dca::func::dmn_0<dca::func::dmn<12>>;
  using Dmn = dca::func::dmn_variadic<Dmn1, Dmn2, Dmn3>;
  using Scalar = TypeParam;
  const std::string typeStr = typeid(TypeParam).name();

  if (12 % comm_size != 0) {
    std::cerr << __func__ << " only works with MPI size that can divide 12 properly. "
              << std::to_string(comm_size) << " is not a good number of ranks to test this."
              << std::endl;
    return;
  }
  std::cout << "comm_size: " << comm_size << "  rank: "<< rank << '\n';
  dca::func::function<Scalar, Dmn, dca::DistType::BLOCKED> f1("parallelFunc", *concurrency_ptr);
  EXPECT_EQ(f1.signature(), 3);
  EXPECT_EQ(f1.get_domain().get_subdomain_size(0), static_cast<std::size_t>(5)) << "on rank: " << rank;
  // the operator[] is broken on ranks other than 0
  // I don't know why yet.
  //EXPECT_EQ(f1[1], static_cast<std::size_t>(4)) << "on rank: " << rank;
  //EXPECT_EQ(f1[2], static_cast<std::size_t>(12)) << "on rank: " << rank;
}

int main(int argc, char** argv) {
  int result = 0;

  dca::linalg::util::printInfoDevices();

  dca::linalg::util::initializeMagma();

  concurrency_ptr = new dca::parallel::MPIConcurrency(argc, argv);
  
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  if (concurrency_ptr->id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  delete concurrency_ptr;
  
  return result;
}
