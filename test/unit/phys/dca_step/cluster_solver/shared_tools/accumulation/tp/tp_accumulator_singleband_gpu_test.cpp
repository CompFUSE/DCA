// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements a no-change test for the two particles accumulation on the GPU.

#include "dca/config/profiler.hpp"

using Scalar = double;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu.hpp"

#include <array>
#include <functional>
#include <string>
#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/distribution/dist_types.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/four_point_type.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

[[maybe_unused]] constexpr bool update_baseline = false;

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_4x4.json";

#ifdef DCA_HAVE_ADIOS2
adios2::ADIOS* adios_ptr;
#endif

#ifdef DCA_HAVE_MPI
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
dca::parallel::MPIConcurrency* concurrency_ptr;
#else
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
dca::parallel::NoConcurrency* concurrency_ptr;
#endif

using ConfigGenerator = dca::testing::AccumulationTest<double>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

using TpAccumulatorGpuSinglebandTest =
    dca::testing::G0Setup<Scalar, dca::testing::LatticeSquare, dca::ClusterSolverId::CT_AUX, input_file>;

unsigned int loop_id = 0;

constexpr bool write_G4s = true;

TEST_F(TpAccumulatorGpuSinglebandTest, Accumulate) {
  dca::linalg::util::initializeMagma();

  const std::array<int, 2> n{27, 24};
  Sample M;
  Configuration config;
  ConfigGenerator::prepareConfiguration(config, M, TpAccumulatorGpuSinglebandTest::BDmn::dmn_size(),
                                        TpAccumulatorGpuSinglebandTest::RDmn::dmn_size(),
                                        parameters_.get_beta(), n);

  using namespace dca::phys;
  std::vector<FourPointType> four_point_channels{FourPointType::PARTICLE_HOLE_TRANSVERSE,
                                                 FourPointType::PARTICLE_HOLE_MAGNETIC,
                                                 FourPointType::PARTICLE_HOLE_CHARGE,
                                                 FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP,
                                                 FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN,
                                                 FourPointType::PARTICLE_PARTICLE_UP_DOWN};

  parameters_.set_four_point_channels(four_point_channels);

  dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::DistType::NONE, dca::linalg::CPU>
      accumulatorHost(data_->G0_k_w_cluster_excluded, parameters_);
  dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::DistType::NONE, dca::linalg::GPU>
      accumulatorDevice(data_->G0_k_w_cluster_excluded, parameters_);
  const int8_t sign = 1;

  accumulatorDevice.resetAccumulation(0);
  accumulatorDevice.accumulate(M, config, sign);
  accumulatorDevice.finalize();

  accumulatorHost.resetAccumulation(0);
  accumulatorHost.accumulate(M, config, sign);
  accumulatorHost.finalize();

#ifdef DCA_HAVE_ADIOS2
  if (write_G4s) {
    dca::io::Writer writer(*adios_ptr, *concurrency_ptr, "ADIOS2", true);
    dca::io::Writer writer_h5(*adios_ptr, *concurrency_ptr, "HDF5", true);

    writer.open_file("tp_single_band_gpu_test_G4.bp");
    writer_h5.open_file("tp_single_band_gpu_test_G4.hdf5");

    parameters_.write(writer);
    parameters_.write(writer_h5);
    data_->write(writer);
    data_->write(writer_h5);

    for (std::size_t channel = 0; channel < accumulatorHost.get_G4().size(); ++channel) {
      std::string channel_str = dca::phys::toString(parameters_.get_four_point_channels()[channel]);
      writer.execute("accumulatorHOST_" + channel_str, accumulatorHost.get_G4()[channel]);
      writer.execute("accumulatorDevice_" + channel_str, accumulatorDevice.get_G4()[channel]);
      writer_h5.execute("accumulatorHOST_" + channel_str, accumulatorHost.get_G4()[channel]);
      writer_h5.execute("accumulatorDevice_" + channel_str, accumulatorDevice.get_G4()[channel]);
    }
    writer.close_file();
    writer_h5.close_file();
  }
#endif
  for (std::size_t channel = 0; channel < accumulatorHost.get_G4().size(); ++channel) {
    const auto diff = dca::func::util::difference(accumulatorHost.get_G4()[channel],
                                                  accumulatorDevice.get_G4()[channel]);
    EXPECT_GT(5e-7, diff.l_inf) << "channel: " << dca::phys::toString(four_point_channels[channel]);
  }
}

int main(int argc, char** argv) {
#ifdef DCA_HAVE_MPI
  dca::parallel::MPIConcurrency concurrency(argc, argv);
  concurrency_ptr = &concurrency;
#else
  dca::parallel::NoConcurrency concurrency(argc, argv);
  concurrency_ptr = &concurrency;
#endif

#ifdef DCA_HAVE_ADIOS2
  // ADIOS expects MPI_COMM pointer or nullptr
  adios2::ADIOS adios("", concurrency_ptr->get(), false);
  adios_ptr = &adios;
#endif

  ::testing::InitGoogleTest(&argc, argv);

  // ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  // delete listeners.Release(listeners.default_result_printer());
  // listeners.Append(new dca::testing::MinimalistPrinter);

  int result = RUN_ALL_TESTS();
  return result;
}
