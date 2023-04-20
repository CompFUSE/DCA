// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter W. Doak     (doakpw@ornl.gov)
//
// This file implements unit tests for two particles accumulation on the GPU.

#include "dca/config/profiler.hpp"
#include "dca/platform/dca_gpu.h"

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
#include <chrono>
#include <thread>
#include "dca/testing/gtest_h_w_warning_blocking.h"
#include "dca/util/integer_division.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/four_point_type.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_4x4_multitransfer_kagome.json";

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

using TpAccumulatorGpuTest =
  dca::testing::G0Setup<Scalar, dca::testing::LatticeKagome, dca::ClusterSolverId::CT_AUX, input_file>;

uint loop_counter = 0;

constexpr bool write_G4s = true;

TEST_F(TpAccumulatorGpuTest, Accumulate) {
  dca::linalg::util::initializeMagma();

  const std::array<int, 2> n{27, 24};
  Sample M;
  Configuration config;
  ConfigGenerator::prepareConfiguration(config, M, TpAccumulatorGpuTest::BDmn::dmn_size(),
                                        TpAccumulatorGpuTest::RDmn::dmn_size(),
                                        parameters_.get_beta(), n);

  using namespace dca::phys;
  std::vector<FourPointType> four_point_channels{
      FourPointType::PARTICLE_HOLE_TRANSVERSE, FourPointType::PARTICLE_HOLE_MAGNETIC,
      FourPointType::PARTICLE_HOLE_CHARGE, FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_UP,
      FourPointType::PARTICLE_HOLE_LONGITUDINAL_UP_DOWN, FourPointType::PARTICLE_PARTICLE_UP_DOWN};
  parameters_.set_four_point_channels(four_point_channels);

  dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::DistType::NONE, dca::linalg::CPU>
      accumulatorHost(data_->G0_k_w_cluster_excluded, parameters_);
  dca::phys::solver::accumulator::TpAccumulator<Parameters, dca::DistType::NONE, dca::linalg::GPU>
      accumulatorDevice(data_->G0_k_w_cluster_excluded, parameters_);
  const int8_t sign = 1;

  accumulatorDevice.resetAccumulation(loop_counter);
  accumulatorDevice.accumulate(M, config, sign);
  accumulatorDevice.finalize();

  accumulatorHost.resetAccumulation(loop_counter);
  accumulatorHost.accumulate(M, config, sign);
  accumulatorHost.finalize();

  ++loop_counter;

#ifdef DCA_HAVE_ADIOS2
  if (write_G4s) {
    dca::io::Writer writer(*adios_ptr, *concurrency_ptr, "ADIOS2", true);
    dca::io::Writer writer_h5(*adios_ptr, *concurrency_ptr, "HDF5", true);

    writer.open_file("tp_gpu_test_G4.bp");
    writer_h5.open_file("tp_gpu_test_G4.hdf5");

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

  std::cout << "blocks: " << dca::util::ceilDiv(int(accumulatorHost.get_G4()[0].size()), 256) << '\n';
  
  for (std::size_t channel = 0; channel < accumulatorHost.get_G4().size(); ++channel) {
    auto diff = dca::func::util::difference(accumulatorHost.get_G4()[channel],
                                            accumulatorDevice.get_G4()[channel]);
    EXPECT_GT(5e-7, diff.l_inf) << "channel: " << dca::phys::toString(four_point_channels[channel]);
  }
}

TEST_F(TpAccumulatorGpuTest, SumToAndFinalize) {
  dca::linalg::util::initializeMagma();

  parameters_.set_four_point_channel(dca::phys::FourPointType::PARTICLE_HOLE_MAGNETIC);

  using Accumulator =
      dca::phys::solver::accumulator::TpAccumulator<G0Setup::Parameters, dca::DistType::NONE,
                                                    dca::linalg::GPU>;
  Accumulator accumulator_sum(data_->G0_k_w_cluster_excluded, parameters_, 0);
  // If the input is not for multiple accumulators this must be manually set or this test can fail
  //accumulator_sum.set_multiple_accumulators(true);
  Accumulator accumulator1(data_->G0_k_w_cluster_excluded, parameters_, 1);
  //accumulator1.set_multiple_accumulators(true);
  Accumulator accumulator2(data_->G0_k_w_cluster_excluded, parameters_, 2);
  //accumulator2.set_multiple_accumulators(true);
  Accumulator accumulator3(data_->G0_k_w_cluster_excluded, parameters_, 3);
  //accumulator3.set_multiple_accumulators(true);


  auto prepare_configuration = [&](auto& M, auto& configuration, const auto& n) {
    ConfigGenerator::prepareConfiguration(M, configuration, TpAccumulatorGpuTest::BDmn::dmn_size(),
                                          TpAccumulatorGpuTest::RDmn::dmn_size(),
                                          parameters_.get_beta(), n);
  };

  const std::array<int, 2> n{3, 4};
  const int8_t sign = -1;
  Sample M1, M2;
  Configuration config1, config2;
  prepare_configuration(config1, M1, n);
  prepare_configuration(config2, M2, n);

  dca::util::OncePerLoopFlag flag;

  int loop_id = loop_counter++;
  std::cout << "loop_id: " << loop_id << '\n';
  accumulator1.resetAccumulation(loop_id, flag);
  accumulator2.resetAccumulation(loop_id, flag);
  accumulator_sum.resetAccumulation(loop_id, flag);

  // This is a bandaid for CI.
  std::chrono::microseconds sleep_time(10000);  
  std::this_thread::sleep_for(sleep_time);

  // there is a data race here between the resets above
  // accumulates because different cuda streams can synchronize the reset and accumulation
  // of the shared G4 data on the GPU.
  accumulator1.accumulate(M1, config1, sign);
  accumulator2.accumulate(M2, config2, sign);
  accumulator1.sumTo(accumulator_sum);
  accumulator2.sumTo(accumulator_sum);  
  
  accumulator_sum.finalize();

  // Reset the G4 on the GPU to zero.

  loop_id = loop_counter++;
  std::cout << "loop_id: " << loop_id << '\n';

  // This is consistent because there is only one accumulator involved and therefore the
  // cuda stream works for synchronization.
  accumulator3.resetAccumulation(loop_id, flag);
  accumulator3.accumulate(M1, config1, sign);
  accumulator3.accumulate(M2, config2, sign);
  accumulator3.finalize();
  
  // auto acc3_it = accumulator3.get_G4()[0].begin();
  // auto acc_sum_it = accumulator_sum.get_G4()[0].begin();
  // auto acc3_end = accumulator3.get_G4()[0].end();

  // int index = 0;
  // while(acc3_it != acc3_end) {
  //   EXPECT_NEAR(acc_sum_it->real(), acc3_it->real(), 1E-4) << "index = " << dca::vectorToString(accumulator3.get_G4()[0].linind_2_subind(index));
  //   ++acc3_it;
  //   ++acc_sum_it;
  //   ++index;
  // }

  const auto diff =
      dca::func::util::difference(accumulator3.get_G4()[0], accumulator_sum.get_G4()[0]);
  EXPECT_GT(5e-7, diff.l_inf);

  if (diff.l_inf < 5e-7) {
    std::cout << "succeed with diff.l_inf = " << diff.l_inf << '\n';
    int index = 10;
    std::cout << "values at " << dca::vectorToString(accumulator3.get_G4()[0].linind_2_subind(index))
	      << " sum: " << accumulator_sum.get_G4()[0](index) << " just_acc: " << accumulator3.get_G4()[0](index) << '\n';
  }
  else {
    int index = 10;
        std::cout << "values at " << dca::vectorToString(accumulator3.get_G4()[0].linind_2_subind(index))
	      << " sum: " << accumulator_sum.get_G4()[0](index) << " just_acc: " << accumulator3.get_G4()[0](index) << '\n';

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
