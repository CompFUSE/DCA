// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This file implements a write read of a largish G4 matrix

#include "dca/platform/dca_gpu.h"

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_cpu.hpp"

#include <array>
#include <map>
#include <string>
#include <vector>

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/io/adios2/adios2_writer.hpp"
#include "dca/io/adios2/adios2_reader.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/four_point_type.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"

using Scalar = double;
#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

int rank, comm_size;
dca::parallel::MPIConcurrency* concurrency_ptr;
adios2::ADIOS* adios_ptr;

constexpr char input_file[] = INPUT_DIR "input_large_G4.json";

using ConfigGenerator = dca::testing::AccumulationTest<double>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

using Scalar = double;

using G4FileIoTest =
  dca::testing::G0Setup<Scalar, dca::testing::LatticeBilayer, dca::ClusterSolverId::CT_AUX, input_file>;

// Since we are not going to put a 1.6G file in the repo this has different logic from tp_accumulator_test.cpp

TEST_F(G4FileIoTest, ReadWrite) {
  dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 0);

  auto fillG4 = [&rng](auto& G4) {
    for (size_t i = 0; i < G4.size(); ++i)
      G4(i) = rng();
  };

  dca::io::ADIOS2Writer writer(concurrency_ptr);
  dca::io::ADIOS2Reader reader(*concurrency_ptr);;

  std::map<dca::phys::FourPointType, std::string> func_names;
  func_names[dca::phys::FourPointType::PARTICLE_HOLE_CHARGE] = "G4_ph_charge";

  const dca::phys::FourPointType g4_channel = dca::phys::FourPointType::PARTICLE_HOLE_CHARGE;

  Data::TpGreensFunction g4_work("G4");
  fillG4(g4_work);

  const std::string self_consistent_large_G4 = "g4_accumulator_test_large_G4.bp";

  writer.open_file(self_consistent_large_G4);
  std::cout << "-- writer.execute name= " << func_names[g4_channel] << std::endl;
  writer.execute(func_names[g4_channel], g4_work);
  writer.close_file();

  Data::TpGreensFunction g4_read("G4");
  reader.open_file(self_consistent_large_G4);
  g4_read.set_name(func_names[g4_channel]);
  reader.execute(g4_read);
  const auto diff = dca::func::util::difference(g4_read, g4_work);
  EXPECT_GT(1e-8, diff.l_inf);
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
    listeners.Append(new dca::testing::
                     MinimalistPrinter);
  }

  adios2::ADIOS adios("", concurrency_ptr->get());
  adios_ptr = &adios;

  result = RUN_ALL_TESTS();

  return result;
}
