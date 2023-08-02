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
#include "dca/util/to_string.hpp"

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_1x1_multi.json";

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

template <typename SCALAR>
struct TpAccumulatorGpuTest : public ::testing::Test {
  using G0Setup = dca::testing::G0SetupBare<SCALAR, dca::testing::LatticeBilayer,
                                            dca::ClusterSolverId::CT_AUX, input_file>;
  virtual void SetUp() {
    host_setup.SetUp();
    gpu_setup.SetUp();
  }

  virtual void TearDown() {}
  G0Setup host_setup;
  G0Setup gpu_setup;
};

uint loop_counter = 0;

constexpr bool write_G4s = true;

using TestTypes = ::testing::Types<double>;
TYPED_TEST_CASE(TpAccumulatorGpuTest, TestTypes);

#define TYPING_PREFACE   using Scalar = TypeParam;\
  using ConfigGenerator = dca::testing::AccumulationTest<Scalar>;\
  using Configuration = typename ConfigGenerator::Configuration;\
  using Sample = typename ConfigGenerator::Sample;


TYPED_TEST(TpAccumulatorGpuTest, Accumulate) {
  TYPING_PREFACE
  
  const std::array<int, 2> n{18, 22};
  Sample M;
  Configuration config;


  ConfigGenerator::prepareConfiguration(config, M,
                                        TpAccumulatorGpuTest<Scalar>::G0Setup::BDmn::dmn_size(),
                                        TpAccumulatorGpuTest<Scalar>::G0Setup::RDmn::dmn_size(),
                                        this->host_setup.parameters_.get_beta(), n);

  using namespace dca::phys;
  std::vector<FourPointType> four_point_channels{
      FourPointType::PARTICLE_HOLE_TRANSVERSE, FourPointType::PARTICLE_HOLE_MAGNETIC,
      FourPointType::PARTICLE_HOLE_CHARGE, FourPointType::PARTICLE_PARTICLE_UP_DOWN};
  this->host_setup.parameters_.set_four_point_channels(four_point_channels);
  this->gpu_setup.parameters_.set_four_point_channels(four_point_channels);

  dca::phys::solver::accumulator::TpAccumulator<decltype(this->host_setup.parameters_),
                                                dca::DistType::NONE, dca::linalg::CPU>
      accumulatorHost(this->host_setup.data_->G0_k_w_cluster_excluded, this->host_setup.parameters_);
  dca::phys::solver::accumulator::TpAccumulator<decltype(this->gpu_setup.parameters_),
                                                dca::DistType::NONE, dca::linalg::GPU>
      accumulatorDevice(this->gpu_setup.data_->G0_k_w_cluster_excluded, this->gpu_setup.parameters_);
  const int8_t sign = 1;

  accumulatorDevice.resetAccumulation(loop_counter);
  accumulatorDevice.accumulate(M, config, sign);
  accumulatorDevice.finalize();

  accumulatorHost.resetAccumulation(loop_counter);
  accumulatorHost.accumulate(M, config, sign);
  accumulatorHost.finalize();

  ++loop_counter;

  // accumulatorDevice.synchronizeStreams();
  
#ifdef DCA_HAVE_ADIOS2
  if (write_G4s) {
    dca::io::Writer writer(*adios_ptr, *concurrency_ptr, "ADIOS2", true);
    dca::io::Writer writer_h5(*adios_ptr, *concurrency_ptr, "HDF5", true);

    writer.open_file("tp_gpu_test_G4.bp");
    writer_h5.open_file("tp_gpu_test_G4.hdf5");

    this->host_setup.parameters_.write(writer);
    this->host_setup.parameters_.write(writer_h5);
    this->host_setup.data_->write(writer);
    this->host_setup.data_->write(writer_h5);

    for (std::size_t channel = 0; channel < accumulatorHost.get_G4().size(); ++channel) {
      std::string channel_str =
          dca::phys::toString(this->host_setup.parameters_.get_four_point_channels()[channel]);
      writer.execute("accumulatorHOST_" + channel_str, accumulatorHost.get_G4()[channel]);
      writer.execute("accumulatorDevice_" + channel_str, accumulatorDevice.get_G4()[channel]);
      writer_h5.execute("accumulatorHOST_" + channel_str, accumulatorHost.get_G4()[channel]);
      writer_h5.execute("accumulatorDevice_" + channel_str, accumulatorDevice.get_G4()[channel]);
    }
    writer_h5.execute("accumulatorDevice_G_0", accumulatorDevice.get_G_Debug()[0]);
    writer_h5.execute("accumulatorDevice_G_1", accumulatorDevice.get_G_Debug()[1]);
    writer_h5.execute("accumulatorHOST_G", accumulatorHost.get_G_Debug());
    writer.close_file();
    writer_h5.close_file();

    // writer_h5.open_file("tp_gpu_test_device_func.hdf5");    
    // this->gpu_setup.parameters_.write(writer_h5);
    // this->gpu_setup.data_->write(writer_h5);
    // writer_h5.close_file();
  }
#endif

  std::cout << "blocks: " << dca::util::ceilDiv(int(accumulatorHost.get_G4()[0].size()), 256) << '\n';

  for (std::size_t channel = 0; channel < accumulatorHost.get_G4().size(); ++channel) {
    auto diff = dca::func::util::difference(accumulatorHost.get_G4()[channel],
                                            accumulatorDevice.get_G4()[channel]);
    EXPECT_GT(5e-7, diff.l_inf) << "channel: " << dca::phys::toString(four_point_channels[channel]);
  }
}

// TYPED_TEST(TpAccumulatorGpuTest, computeM) {
//   TYPING_PREFACE

//   const std::array<int, 2> n{18, 22};
//   Sample M;
//   std::array<dca::linalg::Matrix<Scalar, dca::linalg::GPU>, 2> M_dev;
//   Configuration config;
//   ConfigGenerator::prepareConfiguration(config, M,
//                                         TpAccumulatorGpuTest<Scalar>::G0Setup::BDmn::dmn_size(),
//                                         TpAccumulatorGpuTest<Scalar>::G0Setup::RDmn::dmn_size(),
//                                         this->host_setup.parameters_.get_beta(), n);

//   using namespace dca::phys;
//   std::vector<FourPointType> four_point_channels{
//       FourPointType::PARTICLE_HOLE_TRANSVERSE, FourPointType::PARTICLE_HOLE_MAGNETIC,
//       FourPointType::PARTICLE_HOLE_CHARGE, FourPointType::PARTICLE_PARTICLE_UP_DOWN};
//   this->host_setup.parameters_.set_four_point_channels(four_point_channels);
//   this->gpu_setup.parameters_.set_four_point_channels(four_point_channels);

//   dca::phys::solver::accumulator::TpAccumulator<decltype(this->host_setup.parameters_),
//                                                 dca::DistType::NONE, dca::linalg::CPU>
//       accumulatorHost(this->host_setup.data_->G0_k_w_cluster_excluded, this->host_setup.parameters_);
//   dca::phys::solver::accumulator::TpAccumulator<decltype(this->gpu_setup.parameters_),
//                                                 dca::DistType::NONE, dca::linalg::GPU>
//       accumulatorDevice(this->gpu_setup.data_->G0_k_w_cluster_excluded, this->gpu_setup.parameters_);
//   const int8_t sign = 1;

//   for (int s = 0; s < 2; ++s)
//     M_dev[s].setAsync(M[s], *(accumulatorDevice.get_stream()));

//   accumulatorDevice.resetAccumulation(loop_counter);
//   accumulatorDevice.computeM(M_dev, config);

//   accumulatorHost.resetAccumulation(loop_counter);
//   accumulatorHost.computeM(M, config);

//   Sample M_from_dev;

//   for (int s = 0; s < 2; ++s)
//     M_from_dev[s].setAsync(M_dev[s], *(accumulatorDevice.get_stream()));

//   accumulatorDevice.synchronizeStreams();

//   for (int s = 0; s < 2; ++s)
//   for (int i = 0; i < M[s].nrCols(); ++i)
//     for (int j = 0; j < M[s].nrRows(); ++j) {
//       auto diff = M[s](i, j) - M_from_dev[s](i, j);
//       auto diff_sq = diff * diff;
//       EXPECT_GT(5e-7, diff_sq) << "M[" << s << "](i,j) !=  M_from_dev[" << s << "](i,j) for i:" << i << " j:" << j;
//     }

//   // accumulatorDevice.computeG();
//   // accumulatorHost.computeG();

//   // accumulatorDevice.synchronizeStreams();

//   // using SpGreensFunction = typename decltype(accumulatorHost)::Base::SpGreensFunction;
//   // SpGreensFunction G_from_device("G_from_Device");
//   // dca::util::print_type<SpGreensFunction> print_spg;
//   // print_spg.print(std::cout);

//   // auto leaf_domains = G_from_device.get_domain().get_leaf_domain_sizes();
//   // std::cout << "SpGreensFunction leaf sizes:" << dca::vectorToString(leaf_domains) << '\n';
//   // std::cout << "SpGreensFunction size:" << G_from_device.size() << '\n';
  
//   // using TpComplex = typename decltype(accumulatorHost)::Base::TpComplex;
  
//   // using RMatHost = dca::linalg::ReshapableMatrix<TpComplex, dca::linalg::CPU, dca::config::McOptions::TpAllocator<TpComplex>>;
//   // RMatHost rmat_up = accumulatorDevice.getG(0);
//   // RMatHost rmat_down = accumulatorDevice.getG(1);
//   // std::cout << "rmat nRows: " << rmat_up.nrRows() << " nCols: " << rmat_up.nrCols() << '\n';

//   // using WTpExtDmn = typename decltype(accumulatorHost)::Base::WTpExtDmn;
//   // using KDmn = typename decltype(accumulatorHost)::Base::KDmn;
//   // auto w_size = WTpExtDmn::dmn_size();
//   // auto k_size = KDmn::dmn_size();
//   // for (int w2 = 0; w2 < WTpExtDmn::dmn_size(); ++w2)
//   //   for (int w1 = 0; w1 < WTpExtDmn::dmn_size(); ++w1)
//   //     for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
//   //       for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
//   // 	  auto& g_up = G_from_device(0,0,0,k1,k2,w1,w2);
//   // 	  auto& g_down = G_from_device(0,0,1,k1,k2,w1,w2);
//   // 	  g_up = rmat_up(k1 + k_size * w1, k2 + k_size * w2);
//   // 	  g_down = rmat_down(k1 + k_size * w1, k2 + k_size * w2);
//   // 	}

//   // auto G = accumulatorHost.getG();
  
//   // auto diff = dca::func::util::difference(G_from_device, G);
//   // EXPECT_GT(5e-7, diff.l_inf);
//   // ++loop_counter;
// }

// // TEST_F(TpAccumulatorGpuTest, SumToAndFinalize) {
// //   dca::linalg::util::initializeMagma();

// //   parameters_.set_four_point_channel(dca::phys::FourPointType::PARTICLE_HOLE_MAGNETIC);

// //   using Accumulator =
// //       dca::phys::solver::accumulator::TpAccumulator<G0Setup::Parameters, dca::DistType::NONE,
// //                                                     dca::linalg::GPU>;
// //   Accumulator accumulator_sum(data_->G0_k_w_cluster_excluded, parameters_, 0);
// //   // If the input is not for multiple accumulators this must be manually set or this test can fail
// //   //accumulator_sum.set_multiple_accumulators(true);
// //   Accumulator accumulator1(data_->G0_k_w_cluster_excluded, parameters_, 1);
// //   //accumulator1.set_multiple_accumulators(true);
// //   Accumulator accumulator2(data_->G0_k_w_cluster_excluded, parameters_, 2);
// //   //accumulator2.set_multiple_accumulators(true);
// //   Accumulator accumulator3(data_->G0_k_w_cluster_excluded, parameters_, 3);
// //   //accumulator3.set_multiple_accumulators(true);

// //   auto prepare_configuration = [&](auto& M, auto& configuration, const auto& n) {
// //     ConfigGenerator::prepareConfiguration(M, configuration, TpAccumulatorGpuTest::BDmn::dmn_size(),
// //                                           TpAccumulatorGpuTest::RDmn::dmn_size(),
// //                                           parameters_.get_beta(), n);
// //   };

// //   const std::array<int, 2> n{3, 4};
// //   const int8_t sign = -1;
// //   Sample M1, M2;
// //   Configuration config1, config2;
// //   prepare_configuration(config1, M1, n);
// //   prepare_configuration(config2, M2, n);

// //   dca::util::OncePerLoopFlag flag;

// //   int loop_id = loop_counter++;
// //   std::cout << "loop_id: " << loop_id << '\n';
// //   accumulator1.resetAccumulation(loop_id, flag);
// //   accumulator2.resetAccumulation(loop_id, flag);
// //   accumulator_sum.resetAccumulation(loop_id, flag);

// //   // This is a bandaid for CI.
// //   std::chrono::microseconds sleep_time(10000);
// //   std::this_thread::sleep_for(sleep_time);

// //   // there is a data race here between the resets above
// //   // accumulates because different cuda streams can synchronize the reset and accumulation
// //   // of the shared G4 data on the GPU.
// //   accumulator1.accumulate(M1, config1, sign);
// //   accumulator2.accumulate(M2, config2, sign);
// //   accumulator1.sumTo(accumulator_sum);
// //   accumulator2.sumTo(accumulator_sum);

// //   accumulator_sum.finalize();

// //   // Reset the G4 on the GPU to zero.

// //   loop_id = loop_counter++;
// //   std::cout << "loop_id: " << loop_id << '\n';

// //   // This is consistent because there is only one accumulator involved and therefore the
// //   // cuda stream works for synchronization.
// //   accumulator3.resetAccumulation(loop_id, flag);
// //   accumulator3.accumulate(M1, config1, sign);
// //   accumulator3.accumulate(M2, config2, sign);
// //   accumulator3.finalize();

// //   // auto acc3_it = accumulator3.get_G4()[0].begin();
// //   // auto acc_sum_it = accumulator_sum.get_G4()[0].begin();
// //   // auto acc3_end = accumulator3.get_G4()[0].end();

// //   // int index = 0;
// //   // while(acc3_it != acc3_end) {
// //   //   EXPECT_NEAR(acc_sum_it->real(), acc3_it->real(), 1E-4) << "index = " <<
// //   dca::vectorToString(accumulator3.get_G4()[0].linind_2_subind(index));
// //   //   ++acc3_it;
// //   //   ++acc_sum_it;
// //   //   ++index;
// //   // }

// //   const auto diff =
// //       dca::func::util::difference(accumulator3.get_G4()[0], accumulator_sum.get_G4()[0]);
// //   EXPECT_GT(5e-7, diff.l_inf);

// //   if (diff.l_inf < 5e-7) {
// //     std::cout << "succeed with diff.l_inf = " << diff.l_inf << '\n';
// //     int index = 10;
// //     std::cout << "values at " <<
// //     dca::vectorToString(accumulator3.get_G4()[0].linind_2_subind(index))
// // 	      << " sum: " << accumulator_sum.get_G4()[0](index) << " just_acc: " <<
// // accumulator3.get_G4()[0](index) << '\n';
// //   }
// //   else {
// //     int index = 10;
// //         std::cout << "values at " <<
// //         dca::vectorToString(accumulator3.get_G4()[0].linind_2_subind(index))
// // 	      << " sum: " << accumulator_sum.get_G4()[0](index) << " just_acc: " <<
// // accumulator3.get_G4()[0](index) << '\n';

// //   }
// // }

int main(int argc, char** argv) {
#ifdef DCA_HAVE_MPI
  dca::parallel::MPIConcurrency concurrency(argc, argv);
  concurrency_ptr = &concurrency;
#else
  dca::parallel::NoConcurrency concurrency(argc, argv);
  concurrency_ptr = &concurrency;
#endif

  dca::linalg::util::initializeMagma();

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
