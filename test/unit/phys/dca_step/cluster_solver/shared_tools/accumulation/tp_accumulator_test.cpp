// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements a no-change test for the two point accumulation of a mock configuration.

#include "dca/phys/dca_step/cluster_solver/ctaux/accumulator/tp/accumulator_nonlocal_chi.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/accumulator/tp/accumulator_nonlocal_g.hpp"

#include <array>
#include <string>
#include <vector>
#include <dca/phys/dca_data/dca_data.hpp>
#include "gtest/gtest.h"

#include "dca/function/util/difference.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/four_point_type.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

constexpr bool update_baseline = false;

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/"

constexpr char input_file[] = INPUT_DIR "input_3x3.json";

using ConfigGenerator = dca::testing::AccumulationTest<double>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

using TpAccumulatorTest =
    dca::testing::G0Setup<dca::testing::LatticeBilayer, dca::phys::solver::CT_AUX, input_file>;

TEST_F(TpAccumulatorTest, Accumulate) {
  const std::array<int, 2> n{18, 22};
  Sample M;
  Configuration config;
  ConfigGenerator::prepareConfiguration(config, M, TpAccumulatorTest::BDmn::dmn_size(),
                                        TpAccumulatorTest::RDmn::dmn_size(), parameters_.get_beta(),
                                        n);

  Data::TpGreensFunction G4("G4");
  Data::TpGreensFunction G4_check(G4.get_name());

  dca::phys::solver::ctaux::accumulator_nonlocal_G<Parameters, Data> nonlocal_G_obj(parameters_,
                                                                                    *data_, 0);
  dca::phys::solver::ctaux::accumulator_nonlocal_chi<Parameters, Data> nonlocal_chi_obj(
      parameters_, *data_, 0, G4);

  dca::io::HDF5Writer writer;
  dca::io::HDF5Reader reader;

  const std::string baseline = INPUT_DIR "tp_accumulator_test_baseline.hdf5";

  if (update_baseline)
    writer.open_file(baseline);
  else
    reader.open_file(baseline);

  for (const dca::phys::FourPointType type :
       {dca::phys::PARTICLE_HOLE_TRANSVERSE, dca::phys::PARTICLE_HOLE_MAGNETIC,
        dca::phys::PARTICLE_HOLE_CHARGE, dca::phys::PARTICLE_PARTICLE_UP_DOWN}) {
    G4 = 0;
    parameters_.set_four_point_type(type);

    nonlocal_G_obj.execute(config[0], M[0], config[1], M[1]);
    const int sign = 1;
    nonlocal_chi_obj.execute(sign, nonlocal_G_obj);

    if (update_baseline) {
      G4.set_name("G4_" + toString(type));
      writer.execute(G4);
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
