// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter W. Doak     (doakpw@ornl.gov)
//
// This file implements a no-change test for the two particles accumulation on the GPU with
// the Rashba model.

#include "dca/config/profiler.hpp"
#include <complex>

#include "dca/platform/dca_gpu.h"
#include "dca/util/to_string.hpp"

using Scalar = std::complex<double>;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "dca/phys/dca_step/lattice_mapping/lattice_mapping_sp.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator_gpu.hpp"
#include "dca/function/domains.hpp"

#include <array>
#include <functional>
#include <string>
#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/function/util/difference.hpp"
#include "dca/linalg/util/util_cublas.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/four_point_type.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

constexpr bool write_G4s = true;

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

#define INPUT_DIR \
  DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/"

constexpr char input_file[] = INPUT_DIR "input_222-2_rashba.json";

using ConfigGenerator = dca::testing::AccumulationTest<std::complex<double>>;
using Configuration = ConfigGenerator::Configuration;
using Sample = ConfigGenerator::Sample;

template <typename SCALAR>
struct TpAccumulatorComplexG0GpuTest : public ::testing::Test {
  using G0Setup = dca::testing::G0SetupBare<SCALAR, dca::testing::LatticeRashba,
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

using TestTypes = ::testing::Types<std::complex<double>>;
TYPED_TEST_CASE(TpAccumulatorComplexG0GpuTest, TestTypes);

using namespace dca::phys;

template <class Parameters>
using k_DCA =
    dca::func::dmn_0<domains::cluster_domain<double, Parameters::lattice_type::DIMENSION, domains::CLUSTER,
                                             domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
template <class Parameters>
using k_HOST =
    dca::func::dmn_0<domains::cluster_domain<double, Parameters::lattice_type::DIMENSION, domains::LATTICE_SP,
                                             domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
template <class Parameters, class k_DCA, class k_HOST>
using LatticeMapSpType = latticemapping::lattice_mapping_sp<Parameters, k_DCA, k_HOST>;

#define TYPING_PREFACE                                            \
  using Scalar = TypeParam;                                       \
  using ConfigGenerator = dca::testing::AccumulationTest<Scalar>; \
  using Configuration = typename ConfigGenerator::Configuration;  \
  using Sample = typename ConfigGenerator::Sample;

TYPED_TEST(TpAccumulatorComplexG0GpuTest, Accumulate) {
  TYPING_PREFACE

  const std::array<int, 2> n{23, 23};
  Sample M;
  Configuration config;
  using FourPointType = dca::phys::FourPointType;

  ConfigGenerator::prepareConfiguration(
      config, M, TpAccumulatorComplexG0GpuTest<Scalar>::G0Setup::BDmn::dmn_size(),
      TpAccumulatorComplexG0GpuTest<Scalar>::G0Setup::RDmn::dmn_size(),
      this->host_setup.parameters_.get_beta(), n);
  std::vector<FourPointType> four_point_channels{FourPointType::PARTICLE_PARTICLE_UP_DOWN};

  using namespace dca::phys;
  this->host_setup.parameters_.set_four_point_channels(four_point_channels);
  this->gpu_setup.parameters_.set_four_point_channels(four_point_channels);

  // this->host_setup.data_->initializeSigma("zero");
  // this->gpu_setup.data_->initializeSigma("zero"); //this->gpu_setup.parameters_.get_initial_self_energy());

  using ParametersHost = typename decltype(this->host_setup)::Parameters;
  using ParametersGPU = typename decltype(this->gpu_setup)::Parameters;

  // LatticeMapSpType<ParametersHost,
  // 		   k_DCA<ParametersHost>,
  // 		   k_HOST<ParametersHost>> lattice_mapping_obj_host(this->host_setup.parameters_);
  // auto& host_data = this->host_setup.data_;
  // lattice_mapping_obj_host.execute(host_data->Sigma, host_data->Sigma_lattice_interpolated,
  //                                  host_data->Sigma_lattice_coarsegrained, host_data->Sigma_lattice);

  // LatticeMapSpType<ParametersGPU, k_DCA<ParametersGPU>, k_HOST<ParametersGPU>>
  // lattice_mapping_obj_gpu(this->gpu_setup.parameters_); auto& gpu_data = this->gpu_setup.data_;
  // lattice_mapping_obj_gpu.execute(gpu_data->Sigma, gpu_data->Sigma_lattice_interpolated,
  //                                 gpu_data->Sigma_lattice_coarsegrained, gpu_data->Sigma_lattice);

  dca::phys::solver::accumulator::TpAccumulator<decltype(this->host_setup.parameters_),
                                                dca::DistType::NONE, dca::linalg::CPU>
      accumulatorHost(this->host_setup.data_->G0_k_w_cluster_excluded, this->host_setup.parameters_);
  dca::phys::solver::accumulator::TpAccumulator<decltype(this->gpu_setup.parameters_),
                                                dca::DistType::NONE, dca::linalg::GPU>
      accumulatorDevice(this->gpu_setup.data_->G0_k_w_cluster_excluded, this->gpu_setup.parameters_);
  const std::complex<double> sign = {1.0, 0.0};

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

    writer.open_file("tp_gpu_test_complex_G0_G4.bp");
    writer_h5.open_file("tp_gpu_test_complex_G0_G4.hdf5");

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

#ifndef NDEBUG
    const auto& G_up = accumulatorDevice.get_G_Debug()[0];
    const auto& G_down = accumulatorDevice.get_G_Debug()[1];
    using Parameters = decltype(this->host_setup.parameters_);
    using TpComplex = typename decltype(accumulatorDevice)::TpComplex;
    using HostSpinSepG = dca::linalg::ReshapableMatrix<TpComplex, dca::linalg::CPU,
                                                       dca::linalg::util::PinnedAllocator<TpComplex>>;
    std::array<HostSpinSepG, 2> G_spin_separated{G_up.size(), G_down.size()};

    using WTpExtDmn = dca::func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED>>;
    using KDmn = typename Parameters::KClusterDmn;
    using BDmn = dca::func::dmn_0<domains::electron_band_domain>;
    using SDmn = dca::func::dmn_0<domains::electron_spin_domain>;
    auto& g_all = accumulatorHost.get_G_Debug();

    for (int spin = 0; spin < SDmn::dmn_size(); ++spin) {
      auto& g_this_spin = G_spin_separated[spin];
      auto g_it = g_this_spin.begin();
      for (int w1 = 0; w1 < WTpExtDmn::dmn_size(); ++w1)
        for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1)
          for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1)
            for (int w2 = 0; w2 < WTpExtDmn::dmn_size(); ++w2)
              for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2)
                for (int b2 = 0; b2 < BDmn::dmn_size(); ++b2, ++g_it)
                  *g_it = g_all(b1, b2, spin, k2, k1, w2, w1);
    }

    writer_h5.execute("accumulatorHOST_G_0", G_spin_separated[0]);
    writer_h5.execute("accumulatorHOST_G_1", G_spin_separated[1]);

    for (int i = 0; i < G_up.size().first; ++i)
      for (int j = 0; j < G_up.size().second; ++j) {
        EXPECT_NEAR(G_up(i, j).real(), G_spin_separated[0](i, j).real(), 1E-12)
            << "( " << i << ", " << j << " )";
        EXPECT_NEAR(G_up(i, j).imag(), G_spin_separated[0](i, j).imag(), 1E-12)
            << "( " << i << ", " << j << " )";
        EXPECT_NEAR(G_down(i, j).real(), G_spin_separated[1](i, j).real(), 1E-12)
            << "( " << i << ", " << j << " )";
        EXPECT_NEAR(G_down(i, j).imag(), G_spin_separated[1](i, j).imag(), 1E-12)
            << "( " << i << ", " << j << " )";
      }
#endif

    writer.close_file();
    writer_h5.close_file();
  }
#endif

  std::cout << "blocks: " << dca::util::ceilDiv(int(accumulatorHost.get_G4()[0].size()), 256) << '\n';

  for (std::size_t channel = 0; channel < accumulatorHost.get_G4().size(); ++channel) {
    auto& host_G4 = accumulatorHost.get_G4()[channel];
    auto& device_G4 = accumulatorDevice.get_G4()[channel];
    auto lin_size = host_G4.size();
    int fail_limit = lin_size;
    int fail_count = 0;
    auto& G4_domain_sizes = host_G4.get_domain().get_leaf_domain_sizes();
    for (std::size_t i = 0; i < lin_size; ++i) {
      if (std::abs(host_G4(i).real() - device_G4(i).real()) > 1E-12 ||
          std::abs(host_G4(i).imag() - device_G4(i).imag()) > 1E-12) {
        std::vector<int> fail_index(host_G4.get_domain().get_leaf_domain_sizes().size());
        host_G4.linind_2_subind(i, fail_index);
        ADD_FAILURE() << "host (" << host_G4(i).real() << ", " << host_G4(i).imag()
                      << ") and device (" << device_G4(i).real() << ", " << device_G4(i).imag()
                      << ") G4 fail to match at " << dca::vectorToString(fail_index) << " : "
                      << dca::vectorToString(G4_domain_sizes) << '\n';
        if (++fail_count >= fail_limit)
          break;
      }
      else {
        std::vector<int> success_index(host_G4.get_domain().get_leaf_domain_sizes().size());
        host_G4.linind_2_subind(i, success_index);
      }
    }
    if (fail_count > 0)
      std::cout << "G4 " << fail_count << " of " << lin_size << " elements were incorrect\n";
  }

  for (std::size_t channel = 0; channel < accumulatorHost.get_G4().size(); ++channel) {
    auto diff = dca::func::util::difference(accumulatorHost.get_G4()[channel],
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
