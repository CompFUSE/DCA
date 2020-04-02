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

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator.hpp"

#include <array>
#include <map>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/util/difference.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/four_point_type.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_large_G4.json";

using ConfigGenerator = dca::testing::AccumulationTest<double>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

using G4FileIoTest =
    dca::testing::G0Setup<dca::testing::LatticeBilayer, dca::phys::solver::CT_AUX, input_file>;

// Since we are not going to put a 1.6G file in the repo this has different logic from tp_accumulator_test.cpp

TEST_F(G4FileIoTest, ReadWrite) {
  const std::array<int, 2> n{18, 22};
  
  dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 0);  

  auto fillG4 = [&rng](auto& G4) {
                    for (size_t i = 0; i < G4.size(); ++i)
                      G4(i) = rng();
                };
                      
  dca::io::HDF5Writer writer;
  dca::io::HDF5Reader reader;

  std::map<dca::phys::FourPointType, std::string> func_names;
  func_names[dca::phys::PARTICLE_HOLE_CHARGE] = "G4_ph_charge";

  const dca::phys::FourPointType g4_channel = dca::phys::PARTICLE_HOLE_CHARGE;

  Data::TpGreensFunction g4_work("G4");
  fillG4(g4_work);

  const std::string self_consistent_large_G4 =
      "g4_accumulator_test_large_G4.hdf5";

  writer.open_file(self_consistent_large_G4);
  writer.execute(func_names[g4_channel], g4_work);
  writer.close_file();
  
  Data::TpGreensFunction g4_read("G4");
  reader.open_file(self_consistent_large_G4);
  g4_read.set_name(func_names[g4_channel]);
  reader.execute(g4_read);
  const auto diff = dca::func::util::difference(g4_read, g4_work);
  EXPECT_GT(1e-8, diff.l_inf);
}
