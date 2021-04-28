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
#include "dca/util/to_string.hpp"

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
  using Dmn2 = dca::func::dmn_0<dca::func::dmn<12>>;
  using Dmn3 = dca::func::dmn_0<dca::func::dmn<4>>;
  using Dmn = dca::func::dmn_variadic<Dmn1, Dmn2, Dmn3>;
  using Scalar = TypeParam;
  const std::string typeStr = typeid(TypeParam).name();

  if (4 % comm_size != 0) {
    std::cerr << __func__ << " only works with MPI size that can divide 4 properly. "
              << std::to_string(comm_size) << " is not a good number of ranks to test this."
              << std::endl;
    return;
  }
  std::cout << "comm_size: " << comm_size << "  rank: " << rank << '\n';
  dca::func::function<Scalar, Dmn, dca::DistType::BLOCKED> f1("parallelFunc", *concurrency_ptr);
  EXPECT_EQ(f1.signature(), 3);
  EXPECT_EQ(f1.get_domain().get_subdomain_size(0), static_cast<std::size_t>(5))
      << "on rank: " << rank;
}

TYPED_TEST(BlockDistributedFunctionTest, CheckStartEndSizesOnRanks) {
  int rank = concurrency_ptr->id();
  int comm_size = concurrency_ptr->number_of_processors();

  using Dmn1 = dca::func::dmn_0<dca::func::dmn<5>>;
  using Dmn2 = dca::func::dmn_0<dca::func::dmn<11>>;
  using Dmn3 = dca::func::dmn_0<dca::func::dmn<5>>;
  using Dmn = dca::func::dmn_variadic<Dmn1, Dmn2, Dmn3>;
  using Scalar = TypeParam;
  const std::string typeStr = typeid(TypeParam).name();

  if (4 % comm_size != 0) {
    std::cerr << __func__ << " only works with MPI size that can divide 4 properly. "
              << std::to_string(comm_size) << " is not a good number of ranks to test this."
              << std::endl;
    return;
  }
  std::cout << "comm_size: " << comm_size << "  rank: " << rank << "  ";
  dca::func::function<Scalar, Dmn, dca::DistType::LINEAR> f1("parallelFunc", *concurrency_ptr);
  EXPECT_EQ(f1.signature(), 3);
  EXPECT_EQ(f1.get_domain().get_subdomain_size(0), static_cast<std::size_t>(5))
      << "on rank: " << rank;
  std::cout << dca::vectorToString(f1.get_domain().get_leaf_domain_sizes()) << '\n';

  const std::size_t my_concurrency_id = concurrency_ptr->id();
  const std::size_t my_concurrency_size = concurrency_ptr->number_of_processors();

  std::vector<size_t> sizes(my_concurrency_size, 0);
  std::vector<size_t> starts(my_concurrency_size, 0);
  std::vector<size_t> ends(my_concurrency_size, 0);

  for (int i = 0; i < my_concurrency_size; ++i) {
    size_t local_function_size = dca::util::ceilDiv(f1.get_domain().get_size(), my_concurrency_size);
    size_t residue = f1.get_domain().get_size() % my_concurrency_size;
    size_t start = local_function_size * my_concurrency_id;
    if (my_concurrency_id > (residue - 1)) {
      start -= my_concurrency_id - residue;
      --local_function_size;
    }
    size_t end = start + local_function_size;
    ends[i] = end;
    starts[i] = start;
    if (my_concurrency_id == i)
      sizes[i] = f1.size();
    EXPECT_EQ(start, f1.get_start());
    EXPECT_EQ(end, f1.get_end() + 1);
  }
  concurrency_ptr->sum(sizes);
  size_t total_size = std::accumulate(sizes.begin(), sizes.end(), 0.0);
  EXPECT_EQ(f1.get_domain().get_size(), total_size);
}

TYPED_TEST(BlockDistributedFunctionTest, CheckSplitIndexBlocking) {
  int rank = concurrency_ptr->id();
  size_t comm_size = concurrency_ptr->number_of_processors();

  using Dmn1 = dca::func::dmn_0<dca::func::dmn<5>>;
  using Dmn2 = dca::func::dmn_0<dca::func::dmn<4>>;
  using Dmn3 = dca::func::dmn_0<dca::func::dmn<16>>;
  using Dmn = dca::func::dmn_variadic<Dmn1, Dmn2, Dmn3>;
  using Scalar = TypeParam;
  const std::string typeStr = typeid(TypeParam).name();

  if (4 % comm_size != 0) {
    std::cerr << __func__ << " only works with MPI size that can divide 4 properly. "
              << std::to_string(comm_size) << " is not a good number of ranks to test this."
              << std::endl;
    return;
  }

  std::cout << "comm_size: " << comm_size << "  rank: " << rank << "  ";
  dca::func::function<Scalar, Dmn, dca::DistType::BLOCKED> f1("parallelFunc", *concurrency_ptr);
  EXPECT_EQ(f1.signature(), 3);
  EXPECT_EQ(f1.get_domain().get_subdomain_size(0), static_cast<std::size_t>(5))
      << "on rank: " << rank;
  std::cout << dca::vectorToString(f1.get_domain().get_leaf_domain_sizes()) << '\n';

  std::vector<size_t> sizes(comm_size, 0);

  for (int i = 0; i < comm_size; ++i) {
    size_t local_function_size = dca::util::ceilDiv(f1.get_domain().get_size(), comm_size);
    size_t residue = f1.get_domain().get_size() % comm_size;
    size_t start = local_function_size * rank;
    if (rank > (residue - 1)) {
      start -= rank - residue;
      --local_function_size;
    }
    size_t end = start + local_function_size;
    if (rank == i)
      sizes[i] = f1.size();
    EXPECT_EQ(start, f1.get_start());
    EXPECT_EQ(end, f1.get_end() + 1);
  }
  concurrency_ptr->sum(sizes);
  size_t total_size = std::accumulate(sizes.begin(), sizes.end(), 0.0);
  EXPECT_EQ(f1.get_domain().get_size(), total_size);
  std::cout << "rank(" << rank << ") " << dca::vectorToString(f1.linind_2_subind(f1.get_start()))
            << " : " << dca::vectorToString(f1.linind_2_subind(f1.get_end())) << '\n';
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
