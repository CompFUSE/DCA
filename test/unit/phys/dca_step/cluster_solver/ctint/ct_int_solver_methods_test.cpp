#include "dca/phys/dca_step/cluster_solver/ctint/ctint_cluster_solver.hpp"

#include "gtest/gtest.h"

#include "dca/phys/parameters/parameters.hpp"
#include "test/unit/phys/dca_step/cluster_solver/stub_rng.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/domains/common_domains.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/profiling/null_profiler.hpp"

using RngType = dca::testing::StubRng;
using LatticeType = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using ModelType = dca::phys::models::TightBindingModel<LatticeType>;
using Concurrency = dca::parallel::NoConcurrency;
using Parameters =
    dca::phys::params::Parameters<Concurrency, dca::parallel::NoThreading, dca::profiling::NullProfiler,
                                  ModelType, RngType, dca::phys::solver::CT_INT>;
using Data = dca::phys::DcaData<Parameters>;

struct SolverWrapper
    : public dca::phys::solver::CtintClusterSolver<dca::linalg::CPU, Parameters, false> {
  using BaseClass = dca::phys::solver::CtintClusterSolver<dca::linalg::CPU, Parameters>;

  SolverWrapper(Parameters& parameters_ref, typename BaseClass::Data& data_ref)
      : BaseClass(parameters_ref, data_ref) {}

  using BaseClass::computeSigma;
  using BaseClass::computeG_k_w;
  using Nu = dca::func::dmn_variadic<dca::func::dmn_0<dca::phys::domains::electron_band_domain>,
                                     dca::func::dmn_0<dca::phys::domains::electron_spin_domain>>;
  using Rdmn = typename Parameters::RClusterDmn;
  using Kdmn = typename Parameters::KClusterDmn;
  using Wdmn = typename Parameters::WDmn;

  using F_k_w =
      dca::func::function<std::complex<double>, dca::func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>>;
  void setGandMAnalytic(F_k_w& G, F_k_w& G0, F_k_w& M) const;
};

TEST(CtintClusterSolvertest, ComputeSigma) {
  Concurrency concurrency(1, NULL);
  Parameters parameters("", concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(
      DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctint/inputs/input_single_site.json");
  parameters.update_model();
  parameters.update_domains();

  Data data(parameters);

  SolverWrapper solver(parameters, data);
  SolverWrapper::F_k_w G0, G, M, Sigma;
  solver.setGandMAnalytic(G, G0, M);

  solver.computeSigma(G, G0, Sigma);

  const double U = parameters.get_U();
  const double tollerance = 1e-10;
  for (int s = 0; s < 2; s++)
    for (int w_idx = 0; w_idx < SolverWrapper::Wdmn::dmn_size(); w_idx++) {
      const double w = SolverWrapper::Wdmn::get_elements()[w_idx];
      EXPECT_NEAR(-U * U / (4. * w), Sigma(s, s, 0, w_idx).imag(), tollerance);
      EXPECT_NEAR(0, Sigma(s, s, 0, w_idx).real(), tollerance);
    }
}

void SolverWrapper::setGandMAnalytic(F_k_w& G, F_k_w& G0, F_k_w& M) const {
  const double U = parameters_.get_U();

  for (int i = 0; i < Wdmn::dmn_size(); i++) {
    for (int spin = 0; spin < 2; spin++) {
      const double w = Wdmn::get_elements()[i];
      const std::complex<double> iw(0, w);

      G0(spin, spin, 0, i) = 1. / iw;
      G(spin, spin, 0, i) = 1. / (iw - U * U / (4. * iw));
      M(spin, spin, 0, i) = iw + w * w / (iw - U * U / (4. * iw));
    }
  }
}
