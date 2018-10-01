// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements a no-change test for the two point accumulation on a single band lattice.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator.hpp"

#include <array>
#include "gtest/gtest.h"
#include <string>
#include <vector>

#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

constexpr bool update_baseline = false;

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_3x3.json";

using TpAccumulatorSinglebandTest =
    dca::testing::G0Setup<dca::testing::LatticeSquare, dca::phys::solver::CT_AUX, input_file>;

TEST_F(TpAccumulatorSinglebandTest, Accumulate) {
  using ConfigGenerator = dca::testing::AccumulationTest<Parameters::MC_measurement_scalar_type>;
  using Configuration = ConfigGenerator::Configuration;
  using Sample = ConfigGenerator::Sample;

  const std::array<int, 2> n{17, 17};
  Sample M;
  Configuration config;
  ConfigGenerator::prepareConfiguration(config, M, TpAccumulatorSinglebandTest::BDmn::dmn_size(),
                                        TpAccumulatorSinglebandTest::RDmn::dmn_size(), parameters_.get_beta(),
                                        n);

  Data::TpGreensFunction G4_check("G4");

  dca::io::HDF5Writer writer;
  dca::io::HDF5Reader reader;

  const std::string baseline = INPUT_DIR "tp_accumulator_singleband_test_baseline.hdf5";

  if (update_baseline)
    writer.open_file(baseline);
  else
    reader.open_file(baseline);

  for (const dca::phys::FourPointType type :
       {dca::phys::PARTICLE_HOLE_TRANSVERSE, dca::phys::PARTICLE_HOLE_MAGNETIC,
        dca::phys::PARTICLE_HOLE_CHARGE, dca::phys::PARTICLE_PARTICLE_UP_DOWN}) {
    parameters_.set_four_point_type(type);

    dca::phys::solver::accumulator::TpAccumulator<Parameters> accumulator(
        data_->G0_k_w_cluster_excluded, parameters_);

    const int sign = 1;
    accumulator.accumulate(M, config, sign);
    accumulator.finalize();

    const auto& G4 = accumulator.get_sign_times_G4();

    if (update_baseline) {
      writer.execute("G4_" + toString(type), G4);
    }
    else {
      G4_check.set_name("G4_" + toString(type));
      reader.execute(G4_check);
      const auto diff = dca::func::util::difference(G4, G4_check);
      EXPECT_GT(1e-8, diff.l_inf);
    }
  }

  if (update_baseline)
    writer.close_file();
  else
    reader.close_file();
}
