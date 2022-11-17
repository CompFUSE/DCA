// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// Full DCA computation used to provide input state and reference data to the statistical NiO_test.

#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/math/statistical_testing/function_cut.hpp"
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_ct_hyb_cluster_solver.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/util/git_version.hpp"
#include "dca/distribution/dist_types.hpp"
#include "dca_loop_wrapper.hpp"
#include "NiO_setup.hpp"

constexpr auto SOLVER = dca::ClusterSolverId::CT_INT;

using dca::func::dmn_0;
using dca::func::dmn_variadic;
using RDmn =
    dmn_0<dca::phys::domains::cluster_domain<double, 3, dca::phys::domains::CLUSTER, dca::phys::domains::REAL_SPACE,
                                             dca::phys::domains::BRILLOUIN_ZONE>>;

using SigmaDomain = dca::math::util::SigmaDomain<RDmn>;
using SigmaCutDomain = dmn_variadic<dca::math::util::details::Bdmn, RDmn,
                                    dmn_0<dca::math::util::DmnWcut>, dca::math::util::RealImagDmn>;
using CovarianceDomain = dmn_variadic<SigmaCutDomain, SigmaCutDomain>;
using Concurrency = dca::parallel::MPIConcurrency;
using Model = typename dca::testing::Model;
using Rng = typename dca::testing::RandomNumberGenerator;
using Parameters = typename dca::testing::TestParameters<Concurrency, SOLVER>;
using Data = dca::phys::DcaData<Parameters>;
using ImpuritySolver =
    typename dca::testing::QuantumClusterSolver<SOLVER, Concurrency, dca::linalg::CPU>;
using ClusterSolver = dca::phys::solver::StdThreadQmciClusterSolver<ImpuritySolver>;
using DCACalculation = dca::testing::DcaLoopWrapper<Parameters, Data, ClusterSolver>;

template <class DCACALC, class Parameters, dca::ClusterSolverId SOLVER>
struct MakeOutput;

template <class DCACALC, class Parameters>
struct MakeOutput<DCACALC, Parameters, dca::ClusterSolverId::CT_INT>{
  void operator()(DCACALC& calculation, Parameters& parameters, Concurrency& concurrency) {
    
    const int n_frequencies = 10;
    auto G_k_cut = dca::math::util::bandDiagonal(
        dca::math::util::cutFrequency(calculation.getSolver().local_G_k_w(), n_frequencies));
    auto G_k_Savg(G_k_cut);
    G_k_Savg.set_name("G_k_w");
    concurrency.sum_and_average(G_k_Savg);
  }
};

template <class DCACALC, class Parameters>
struct MakeOutput<DCACALC, Parameters, dca::ClusterSolverId::SS_CT_HYB>  {
  void operator()(DCACALC& calculation, Parameters& parameters, Concurrency& concurrency) {
    const int n_frequencies = 10;
    auto GScut = dca::math::util::bandDiagonal(
        dca::math::util::cutFrequency(calculation.getSolver().local_GS_r_w(), n_frequencies));
    auto GSavg(GScut);
    GSavg.set_name("GS_r_w");
    concurrency.sum_and_average(GSavg);

    dca::func::function<double, CovarianceDomain> cov("Covariance");
    concurrency.computeCovariance(GScut, GSavg, cov);
    // write the result
    if (concurrency.id() == 0) {
      dca::io::HDF5Writer writer;
      writer.open_file("NiO_covariance_output.hdf5");
      writer.open_group("functions");
      writer.execute(GSavg);
      writer.execute(cov);
      writer.close_group();
      writer.open_group("parameters");
      auto mpi_size = concurrency.number_of_processors();
      writer.execute("measurments_per_rank", parameters.get_measurements().back() / mpi_size);
      writer.execute("nodes", concurrency.number_of_processors());
      writer.close_group();
      writer.close_file();
    }
  }
};

int main(int argc, char** argv) {
  const std::string input_file(TEST_DIRECTORY "input_NiO_ctint.json");

  Concurrency concurrency(argc, argv);
  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
  // override file input for file names
  parameters.set_t_ij_file_name(TEST_DIRECTORY "t_ij_NiO.txt");
  parameters.set_U_ij_file_name(TEST_DIRECTORY "U_ij_NiO_8_lit.txt");

  parameters.update_model();
  parameters.update_domains();
  Data dca_data(parameters);
  dca_data.initialize();
  const int id = concurrency.id();

  DCACalculation calculation(parameters, dca_data, concurrency);

  calculation.initialize();
  calculation.execute();
  calculation.finalize();

  // do one last loop stopping before MC integration and write the state
  calculation.performPreIntegrationSteps();
  if (id == concurrency.last()) {
    calculation.write();
  }
  // do one more iteration and compute covariance
  calculation.performIntegrationStep();
  MakeOutput<decltype(calculation), decltype(parameters), SOLVER> output;
  output(calculation, parameters, concurrency);
}
