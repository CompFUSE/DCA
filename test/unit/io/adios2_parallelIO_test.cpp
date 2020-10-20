// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides specific tests for the ADIOS2 reader and writer.

#include "dca/io/adios2/adios2_reader.hpp"
#include "dca/io/adios2/adios2_writer.hpp"
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
#include "dca/parallel/util/get_workload.hpp"
#include "dca/testing/minimalist_printer.hpp"

#include <array>
#include <complex>
#include <string>
#include <vector>
#include <typeinfo>  // operator typeid
#include <iostream>
#include <sstream>

#include "gtest/gtest.h"

int rank, comm_size;
dca::parallel::MPIConcurrency* concurrency_ptr;

template <typename Scalar>
class ADIOS2ParallelIOTest : public ::testing::Test {};
// using TestTypes = ::testing::Types<float, std::complex<double>>;
using TestTypes = ::testing::Types<float>;

TYPED_TEST_CASE(ADIOS2ParallelIOTest, TestTypes);

template <class T>
std::string VectorToString(const std::vector<T>& v) {
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < v.size(); ++i) {
    if (i != 0)
      ss << ",";
    ss << v[i];
  }
  ss << "]";
  return ss.str();
}

TYPED_TEST(ADIOS2ParallelIOTest, FunctionReadWrite) {
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

  dca::func::function<Scalar, Dmn> f1("parallelFunc", *concurrency_ptr);
  size_t dmn_size = 1;
  for (int l = 0; l < f1.signature(); ++l)
    dmn_size *= f1[l];

  int val = rank * dmn_size;

  uint64_t start = 0;
  uint64_t end = 0;

  
  
  std::cout << "Rank:" << rank << "  start: " << start << "  end:" << end << '\n';
  // only set this ranks values
  for (int i = 0; i <= end - start; ++i)
    f1.data()[i] = ++val;

  // get the N-dimensional decomposition
  std::vector<int> subind_start = f1.linind_2_subind(start);
  std::vector<int> subind_end = f1.linind_2_subind(end);
  std::cout << "rank: " << rank << " start lin = " << start
            << " subindicies = " << VectorToString(subind_start) << " end lin = " << end
            << " subindicies = " << VectorToString(subind_end) << std::endl;
  EXPECT_EQ(start, static_cast<uint64_t>(rank) * (dmn_size / comm_size));
  EXPECT_EQ(end - start + 1, dmn_size / comm_size);

  const std::string fname("ADIOS2ParallelIOTest_" + typeStr + ".bp");
  {
    dca::io::ADIOS2Writer writer(concurrency_ptr, "", true);
    writer.open_file(fname, true);

    // Because the caller needs to know if its function is distributed or not we will assume this is
    // so for the API as well. in the future I think something more sophisticated needs to be done
    // and the function will need to know its distribution, but for now we distribute only over the
    // fastest index

    writer.execute(f1, subind_start, subind_end);

    writer.close_file();
  }
  {
    // Read test file.
    if (!rank) {
      std::cout << " Read back data with 3D selection " << std::endl;
    }

    dca::io::ADIOS2Reader reader(concurrency_ptr, "", true);
    reader.open_file(fname);

    dca::func::function<Scalar, Dmn> f2("parallelFunc");

    EXPECT_TRUE(reader.execute(f2, subind_start, subind_end));

    /* TODO: This should be working on every rank */
    if (!rank) {
      for (int i = start; i < end; ++i)
        EXPECT_EQ(f1(i), f2(i));
    }

    if (!rank) {
      std::cout << " Read back data with linear 1D selection " << std::endl;
    }

    /* We can use the linear index to read back the 3D array */
    EXPECT_TRUE(reader.execute(f2, start, end));
    /* TODO: This should be working on every rank */
    if (!rank) {
      for (int i = start; i < end; ++i)
        EXPECT_EQ(f1(i), f2(i));
    }

    reader.close_file();
  }
}

TYPED_TEST(ADIOS2ParallelIOTest, FunctionReadWriteLinear) {
  using Dmn1 = dca::func::dmn_0<dca::func::dmn<5>>;
  using Dmn2 = dca::func::dmn_0<dca::func::dmn<4>>;
  using Dmn3 = dca::func::dmn_0<dca::func::dmn<2>>;
  using Dmn = dca::func::dmn_variadic<Dmn1, Dmn2, Dmn3>;
  using Scalar = TypeParam;
  const std::string typeStr = typeid(TypeParam).name();

  dca::func::function<Scalar, Dmn> f1("parallelFunc", *concurrency_ptr);
  size_t dmn_size = 1;
  for (int l = 0; l < f1.signature(); ++l)
    dmn_size *= f1[l];

  int val = rank * dmn_size;

  uint64_t start = 0;
  uint64_t end = 0;
  // This returns the linearized bounds of the function for a rank.
  dca::parallel::util::getComputeRange(concurrency_ptr->id(), concurrency_ptr->number_of_processors(),
                                       static_cast<uint64_t>(f1.size()), start, end);

  // only set this ranks values
  for (int i = start; i <= end; ++i)
    f1.data()[i] = ++val;

  std::cout << "rank: " << rank << " start lin = " << start << " end lin = " << end << std::endl;

  const std::string fname("ADIOS2ParallelIOTest_linear_" + typeStr + ".bp");
  {
    dca::io::ADIOS2Writer writer(concurrency_ptr, "", true);
    writer.open_file(fname, true);

    // Because the caller needs to know if its function is distributed or not we will assume this is
    // so for the API as well. in the future I think something more sophisticated needs to be done
    // and the function will need to know its distribution, but for now we distribute only over the
    // fastest index

    writer.execute(f1, start, end);

    writer.close_file();
  }
  {
    // Read test file.
    if (!rank) {
      std::cout << " Read back data with linear 1D selection " << std::endl;
    }
    dca::io::ADIOS2Reader reader(concurrency_ptr, "", true);
    reader.open_file(fname);

    dca::func::function<Scalar, Dmn> f2("parallelFunc");

    EXPECT_TRUE(reader.execute(f2, start, end));

    /* TODO: This should be working on every rank */
    if (!rank) {
      for (int i = start; i < end; ++i)
        EXPECT_EQ(f1(i), f2(i));
    }

    if (!rank) {
      std::cout << " Attempt to read back data with 3D selection should fail " << std::endl;
    }

    // get the N-dimensional decomposition, which is entirely wrong for this test
    std::vector<int> subind_start = f1.linind_2_subind(start);
    std::vector<int> subind_end = f1.linind_2_subind(end);

    /* 1D array on disk cannot be read back with 3D selection */
    EXPECT_FALSE(reader.execute(f2, subind_start, subind_end));

    reader.close_file();
  }
}

int main(int argc, char** argv) {
  int result = 0;

  dca::parallel::MPIConcurrency concurrency(argc, argv);
  rank = concurrency.id();
  comm_size = concurrency.number_of_processors();
  concurrency_ptr = &concurrency;

  ::testing::InitGoogleTest(&argc, argv);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  if (rank != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
