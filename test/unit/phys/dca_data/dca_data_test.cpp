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

#include <functional>
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
          const char* input_name = dca_data_test_input, DistType DT = DistType::NONE,
          class NUMTRAITS = dca::NumericalTraits<dca::util::RealAlias<Scalar>, Scalar>>
struct DCADataSetup {
  using LatticeType = Lattice;
  using Model = phys::models::TightBindingModel<Lattice>;
  using RngType = testing::StubRng;
  using Parameters =
      phys::params::Parameters<Concurrency, parallel::NoThreading, profiling::NullProfiler, Model,
                               RngType, solver_name, NUMTRAITS>;
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

constexpr char dca_data_test_adios2_input[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_data/input_adios2.json";

using DCADataTest = DCADataTestWrapper<DCADSetup<read_sigma_input_file>>;

TEST_F(DCADataTest, ReadSigma) {
  dca_setup_.data_->initialize();
  std::cout << "initial_self_energy:" << dca_setup_.parameters_->get_initial_self_energy();
  dca_setup_.data_->initializeSigma("dca_data_test.hdf5");
  EXPECT_EQ(dca_setup_.parameters_->get_chemical_potential(), 3.0);
  decltype(dca_setup_.data_->Sigma) check_sigma("check_sigma");
  check_sigma = 3.0;
  EXPECT_EQ(dca_setup_.data_->Sigma, check_sigma);
}

TEST_F(DCADataTest, ReadSigmaLegacy) {
  dca_setup_.data_->initialize();
  dca_setup_.data_->initializeSigma("dca_data_test_legacy.hdf5");
  EXPECT_EQ(dca_setup_.parameters_->get_chemical_potential(), 4.0);
  decltype(dca_setup_.data_->Sigma) check_sigma("check_sigma");
  check_sigma = 4.0;
  EXPECT_EQ(dca_setup_.data_->Sigma, check_sigma);
}

#ifdef DCA_HAVE_ADIOS2
TEST_F(DCADataTest, ReadSigmaAdios2) {
  dca_setup_.data_->initialize();
  dca_setup_.data_->initializeSigma("dca_data_test.bp");
  EXPECT_EQ(dca_setup_.parameters_->get_chemical_potential(), 3.0);
  decltype(dca_setup_.data_->Sigma) check_sigma("check_sigma");
  check_sigma = 3.0;
  EXPECT_EQ(dca_setup_.data_->Sigma, check_sigma);
}

TEST_F(DCADataTest, ReadSigmaAdios2ExtraPot) {
  dca_setup_.data_->initialize();
  dca_setup_.data_->initializeSigma("dca_data_test_extra_cpot.bp");
  EXPECT_EQ(dca_setup_.parameters_->get_chemical_potential(), 3.0);
  decltype(dca_setup_.data_->Sigma) check_sigma("check_sigma");
  check_sigma = 3.0;
  EXPECT_EQ(dca_setup_.data_->Sigma, check_sigma);
}

TEST_F(DCADataTest, ReadSigmaAdios2AnotherPot) {
  dca_setup_.data_->initialize();
  dca_setup_.data_->initializeSigma("dca_data_test_another_cpot.bp");
  EXPECT_EQ(dca_setup_.parameters_->get_chemical_potential(), 4.0);
  decltype(dca_setup_.data_->Sigma) check_sigma("check_sigma");
  check_sigma = 4.0;
  EXPECT_EQ(dca_setup_.data_->Sigma, check_sigma);
}

#endif

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
  adios2::ADIOS adios("", concurrency_ptr->get());
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
    dca_data.data_->Sigma = 1.0;
    std::string format{"HDF5"};
    dca::io::Writer<Concurrency> writer(*concurrency_ptr, format);
    writer.open_file("dca_data_test.hdf5", true);
    writer.begin_step();
    dca_data.data_->write(writer);
    dca_data.loop_data_->write(writer, *concurrency_ptr);
    writer.end_step();
    dca_data.loop_data_->chemical_potential(1) = 2.0;
    dca_data.loop_data_->last_completed_iteration = 1;
    dca_data.data_->Sigma = 2.0;
    writer.begin_step();
    dca_data.data_->write(writer);
    dca_data.loop_data_->write(writer, *concurrency_ptr);
    writer.end_step();
    dca_data.loop_data_->chemical_potential(2) = 3.0;
    dca_data.loop_data_->last_completed_iteration = 2;
    dca_data.data_->Sigma = 3.0;
    writer.begin_step();
    dca_data.data_->write(writer);
    dca_data.loop_data_->write(writer, *concurrency_ptr);
    writer.end_step();
    writer.close_file();

    // Fake legacy hdf5 file by omitting steps
    dca::io::Writer<Concurrency> writer_legacy(*concurrency_ptr, format);
    writer_legacy.open_file("dca_data_test_legacy.hdf5", true);
    dca_data.loop_data_->chemical_potential(3) = 4.0;
    dca_data.loop_data_->last_completed_iteration = 3;
    dca_data.data_->Sigma = 4.0;
    dca_data.data_->write(writer_legacy);
    dca_data.loop_data_->write(writer_legacy, *concurrency_ptr);
    auto& hdf5_writer = std::get<dca::io::HDF5Writer>(writer_legacy.getUnderlying());
    // this writes without a top level steps variable
    hdf5_writer.legacy_close_file();
  }
#ifdef DCA_HAVE_ADIOS2

  {
    // Make adioss file for testing.
    DCADSetup<dca_data_test_adios2_input> dca_data;
    dca_data.SetUp(concurrency_ptr);
    dca_data.data_->initialize();
    dca_data.loop_data_->chemical_potential(0) = 1.0;
    dca_data.loop_data_->last_completed_iteration = 0;
    dca_data.data_->Sigma = 1.0;
    std::string format{"ADIOS2"};
    dca::io::Writer<Concurrency> writer(*concurrency_ptr, format);
    writer.open_file("dca_data_test.bp", true);
    // insure chemical potentials past the last complete iteration are ignored.
    dca::io::Writer<Concurrency> writer_extra_cpot(*concurrency_ptr, format);
    writer_extra_cpot.open_file("dca_data_test_extra_cpot.bp", true);
    // insure chemical potentials past the last complete iteration are ignored.
    dca::io::Writer<Concurrency> writer_another_cpot(*concurrency_ptr, format);
    writer_another_cpot.open_file("dca_data_test_another_cpot.bp", true);

    std::vector<std::reference_wrapper<dca::io::Writer<Concurrency>>> writers{
        writer, writer_extra_cpot, writer_another_cpot};

    auto begin_steps = [](dca::io::Writer<Concurrency>& writer) { writer.begin_step(); };
    auto writes = [&dca_data](dca::io::Writer<Concurrency>& writer) {
      dca_data.data_->write(writer);
      dca_data.loop_data_->write(writer, *concurrency_ptr);
    };
    auto end_steps = [](dca::io::Writer<Concurrency>& writer) { writer.end_step(); };

    std::for_each(writers.begin(), writers.end(), begin_steps);
    std::for_each(writers.begin(), writers.end(), writes);
    std::for_each(writers.begin(), writers.end(), end_steps);

    dca_data.loop_data_->chemical_potential(1) = 2.0;
    dca_data.loop_data_->last_completed_iteration = 1;
    dca_data.data_->Sigma = 2.0;
    std::for_each(writers.begin(), writers.end(), begin_steps);
    std::for_each(writers.begin(), writers.end(), writes);
    std::for_each(writers.begin(), writers.end(), end_steps);

    dca_data.loop_data_->chemical_potential(2) = 3.0;
    dca_data.loop_data_->last_completed_iteration = 2;
    dca_data.data_->Sigma = 3.0;
    std::for_each(writers.begin(), writers.end(), begin_steps);
    std::for_each(writers.begin(), writers.end(), writes);
    std::for_each(writers.begin(), writers.end(), end_steps);

    dca_data.loop_data_->chemical_potential(3) = 4.0;
    dca_data.data_->Sigma = 4.0;
    writer_extra_cpot.begin_step();
    dca_data.data_->write(writer_extra_cpot);
    dca_data.loop_data_->write(writer_extra_cpot, *concurrency_ptr);
    writer_extra_cpot.end_step();
    writer_extra_cpot.close_file();

    dca_data.loop_data_->last_completed_iteration = 3;
    writer_another_cpot.begin_step();
    dca_data.data_->write(writer_another_cpot);
    dca_data.loop_data_->write(writer_another_cpot, *concurrency_ptr);
    writer_another_cpot.end_step();
    writer_another_cpot.close_file();
  }
#endif

  int result = RUN_ALL_TESTS();
  return result;
}
