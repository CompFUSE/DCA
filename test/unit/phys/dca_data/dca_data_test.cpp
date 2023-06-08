// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file provides unit tests for dca_data, thus far these just cover specific cases where there
// have been regressions

#include "dca/testing/gtest_h_w_warning_blocking.h"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "test/unit/phys/dca_step/cluster_solver/stub_rng.hpp"
#include "dca/profiling/null_profiler.hpp"

#ifdef DCA_HAVE_ADIOS2
adios2::ADIOS* adios_ptr;
#endif

#ifdef DCA_HAVE_MPI
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
using Concurrency = dca::parallel::MPIConcurrency;
#else
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
using Concurrency = dca::parallel::NoConcurrency;
#endif
Concurrency* concurrency_ptr;

namespace dca {
namespace testing {
// dca::testing::

using LatticeSquare = phys::models::square_lattice<phys::domains::D4>;

constexpr char dca_data_test_input[] = DCA_SOURCE_DIR "/test/unit/phys/dca_data/input_hdf5.json";

template <typename Scalar, class Concurrency, class Lattice = LatticeSquare,
          ClusterSolverId solver_name = ClusterSolverId::CT_AUX,
          const char* input_name = dca_data_test_input, DistType DT = DistType::NONE>
struct DCADataSetup {
  using LatticeType = Lattice;
  using Model = phys::models::TightBindingModel<Lattice>;
  using RngType = testing::StubRng;
  using Parameters = phys::params::Parameters<Concurrency, parallel::NoThreading,
                                              profiling::NullProfiler, Model, RngType, solver_name>;
  using Data = phys::DcaData<Parameters, DT>;
  using LoopData = phys::DcaLoopData<Parameters>;

  Concurrency* concurrency_;
  std::unique_ptr<Parameters> parameters_;
  std::unique_ptr<Data> data_;
  std::unique_ptr<LoopData> loop_data_;

  void SetUp(Concurrency* concurrency) {
    concurrency_ = concurrency;
    parameters_ = std::make_unique<Parameters>("", *concurrency);
    try {
      parameters_->template read_input_and_broadcast<io::JSONReader>(input_name);
    }
    catch (const std::exception& r_w) {
      throw std::runtime_error(r_w.what());
    }
    catch (...) {
      throw std::runtime_error("Input parsing failed!");
    }
    parameters_->update_model();
    static bool domain_initialized = false;
    if (!domain_initialized) {
      parameters_->update_domains();
      domain_initialized = true;
    }
    data_ = std::make_unique<Data>(*parameters_);
    loop_data_ = std::make_unique<LoopData>();
  }

  void TearDown() {}
};

}  // namespace testing
}  // namespace dca

template <const char* INPUTFILE = dca::testing::dca_data_test_input>
using DCADSetup = dca::testing::DCADataSetup<double, Concurrency, dca::testing::LatticeSquare,
                                             dca::ClusterSolverId::CT_AUX, INPUTFILE>;

template <class DCADATASETUP>
struct DCADataTestWrapper : public ::testing::Test {
  DCADATASETUP dca_setup_;
  virtual void SetUp() {
    dca_setup_.SetUp(concurrency_ptr);
  }
  virtual void TearDown() {
    dca_setup_.TearDown();
  }
};

constexpr char read_sigma_input_file[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_data/input_read_hdf5.json";

using DCADataTest = DCADataTestWrapper<DCADSetup<read_sigma_input_file>>;

TEST_F(DCADataTest, ReadSigma) {
  dca_setup_.data_->initialize();
  std::cout << "initial_self_energy:" << dca_setup_.parameters_->get_initial_self_energy();
  dca_setup_.data_->initializeSigma("dca_data_test.hdf5");
  EXPECT_EQ(dca_setup_.parameters_->get_chemical_potential(), 3.0);
}

TEST_F(DCADataTest, ReadSigmaLegacy) {
  dca_setup_.data_->initialize();
  dca_setup_.data_->initializeSigma("dca_data_test_legacy.hdf5");
  EXPECT_EQ(dca_setup_.parameters_->get_chemical_potential(), 4.0);
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

  {
    // Make hdf file for testing.
    DCADSetup<> dca_data;
    dca_data.SetUp(concurrency_ptr);
    dca_data.data_->initialize();
    dca_data.loop_data_->chemical_potential(0) = 1.0;
    dca_data.loop_data_->last_completed_iteration = 0;
    std::string format{"HDF5"};
    dca::io::Writer<Concurrency> writer(
#ifdef DCA_HAVE_ADIOS2
        *adios_ptr,
#endif
        *concurrency_ptr, format);
    writer.open_file("dca_data_test.hdf5", true);
    writer.begin_step();
    dca_data.data_->write(writer);
    dca_data.loop_data_->write(writer, *concurrency_ptr);
    writer.end_step();
    dca_data.loop_data_->chemical_potential(1) = 2.0;
    dca_data.loop_data_->last_completed_iteration = 1;
    writer.begin_step();
    dca_data.data_->write(writer);
    dca_data.loop_data_->write(writer, *concurrency_ptr);
    writer.end_step();
    dca_data.loop_data_->chemical_potential(2) = 3.0;
    dca_data.loop_data_->last_completed_iteration = 2;
    writer.begin_step();
    dca_data.data_->write(writer);
    dca_data.loop_data_->write(writer, *concurrency_ptr);
    writer.end_step();
    writer.close_file();

	// Fake legacy hdf5 file by omitting steps
    dca::io::Writer<Concurrency> writer_legacy(
#ifdef DCA_HAVE_ADIOS2
        *adios_ptr,
#endif
        *concurrency_ptr, format);
    writer_legacy.open_file("dca_data_test_legacy.hdf5", true);
    dca_data.loop_data_->chemical_potential(3) = 4.0;
    dca_data.loop_data_->last_completed_iteration = 3;
    dca_data.data_->write(writer_legacy);
    dca_data.loop_data_->write(writer_legacy, *concurrency_ptr);
    auto& hdf5_writer = std::get<dca::io::HDF5Writer>(writer_legacy.getUnderlying());
    // this writes without a top level steps variable
    hdf5_writer.legacy_close_file();
  }
  int result = RUN_ALL_TESTS();
  return result;
}
