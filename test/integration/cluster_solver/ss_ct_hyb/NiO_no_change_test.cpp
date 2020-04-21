// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// No change test for the CT-HYB solver using  a NiO lattice.

#include "gtest/gtest.h"

#include "dca/config/threading.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_cluster_solver.hpp"
#include "dca/phys/models/material_hamiltonians/material_lattice.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/util/git_version.hpp"

constexpr int update_baseline = false;

const std::string test_directory = DCA_SOURCE_DIR "/test/integration/cluster_solver/ss_ct_hyb/";

TEST(Ni0NoChangeTest, GreensFunction) {
  using Concurrency = dca::parallel::MPIConcurrency;
  Concurrency concurrency(0, nullptr);

  const int id = concurrency.id();

  if (id == 0)
    dca::util::GitVersion::print();

  using Model = dca::phys::models::TightBindingModel<dca::phys::models::material_lattice<
      dca::phys::models::NiO_unsymmetric, dca::phys::domains::no_symmetry<3>>>;
  using Rng = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
  using TestParameters =
      dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler,
                                    Model, Rng, dca::phys::solver::SS_CT_HYB>;
  using Data = dca::phys::DcaData<TestParameters>;
  using ImpuritySolver = dca::phys::solver::StdThreadQmciClusterSolver<
      dca::phys::solver::SsCtHybClusterSolver<dca::linalg::CPU, TestParameters, Data>>;

  TestParameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(test_directory + "input_NiO.json");
  // override file input for file names
  parameters.set_t_ij_file_name(test_directory + "t_ij_NiO.txt");
  parameters.set_U_ij_file_name(test_directory + "U_ij_NiO_8_lit.txt");
  parameters.update_model();
  parameters.update_domains();

  Data data(parameters);
  // initialize H only. G0 is read from file afterwards.
  data.initializeH0_and_H_i();

  // Read and broadcast the rest of the initialization from full DCA results.
  if (id == 0) {
    dca::io::HDF5Reader reader;
    reader.open_file(test_directory + "NiO_coarse_grained.hdf5");
    reader.open_group("functions");
    reader.execute(data.G0_k_w);
    reader.execute(data.G0_r_t);
    reader.execute(data.G_k_w);
    reader.close_group();
    reader.open_group("additional_functions");
    reader.execute(data.Sigma_cluster);
    reader.close_group();
    reader.close_file();
  }
  concurrency.broadcast(data.G0_k_w);
  concurrency.broadcast(data.G0_r_t);
  concurrency.broadcast(data.G_k_w);
  concurrency.broadcast(data.Sigma_cluster);

  ImpuritySolver solver(parameters, data);
  solver.initialize(0);
  solver.integrate();

  dca::phys::DcaLoopData<TestParameters> loop_data;
  solver.finalize(loop_data);

  const std::string baseline_filename = test_directory + "NiO_baseline.hdf5";

  if (id == 0) {
    if (!update_baseline) {
      dca::io::HDF5Reader reader;
      reader.open_file(baseline_filename);
      reader.open_group("functions");

      Data::SpGreensFunction G_expected(data.G_k_w.get_name());
      reader.execute(G_expected);

      auto err = dca::func::util::difference(data.G_k_w, G_expected);
      EXPECT_LE(err.l_inf, 1e-7);
    }
    else {
      std::cout << "\nProcessor " << id << " is writing data " << std::endl;
      dca::io::HDF5Writer writer;
      writer.open_file(baseline_filename);
      writer.open_group("functions");
      writer.execute(data.G_k_w);
      writer.close_file();
    }
  }
}
