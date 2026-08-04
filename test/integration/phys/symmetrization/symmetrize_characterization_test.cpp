// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Tyler Sax (tylersax@gmail.com)
//
// Characterization harness for single-particle cluster symmetrization.
//
// This is the regression net for Hamiltonian-derived-symmetrization. It
// generalizes symmetrize_test.cpp from one trivial-symmetry case (twoband_chain x
// no_symmetry<2>, an identity-only group that makes every cluster loop a structural
// no-op) into a per-(model, operation) characterization across lattices paired with
// their real point groups.
//
// Oracle: G0 = [iw + mu - H0]^-1 is deterministic and exactly invariant under every
// symmetry S of H0, so symmetrizing G0 must be a no-op. Four tests exercise this:
//   1. DomainActivation: assert the cluster realized the expected number of symmetry
//      ops, so the no-op checks below have a non-trivial group to exercise.
//   2. CoarseNoOp: run the real production Symmetrize::execute on a copy of G0 and diff
//      against the original at ~1e-12, over all four representations (k/r x w/t). This
//      is the authoritative regression net -- it runs the whole production pipeline,
//      which is where the multi-orbital bugs live.
//   3 & 4. PerOpMapMomentum / PerOpMapRealSpace: a test-local applicator that, for one
//      operation index, applies today's cluster-symmetry table relation to G0 and checks
//      invariance. They mirror executeCluster, so they characterize today's
//      representation+code per operation -- isolating which op's encoding is to blame.
//
// Limitation: A no-op on already-symmetric G0 detects operations whose encoding
// is wrong, but does not detect symmetry operations that are silently missing.
//
// Test architecture: One binary per case. The cluster_symmetry<cluster_domain<...>> table
// and all DCA domains are singletons keyed on (scalar, D, role-NAME, rep, shape) -- not on
// the lattice -- so two lattices of the same dimension alias the same singletons and a
// second update_domains() would clobber the first. Each case is therefore built into its
// own binary, selected by THIS_TEST_DEFINES via CMake TEST_DEFINES (same idiom as
// tp_accumulator_gpu_test.cpp).

#include "dca/platform/dca_gpu.h"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/dca_step/symmetrization/derive_point_group.hpp"
#include "dca/phys/dca_step/symmetrization/solve_orbital_op_signs.hpp"

#include <algorithm>
#include <complex>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/config/threading.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/function/util/real_complex_conversion.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/dca_algorithms/compute_free_greens_function.hpp"
#include "dca/phys/domains/cluster/cluster_symmetry.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/holohedries_2d.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/phys/models/analytic_hamiltonians/Kagome_hubbard.hpp"
#include "dca/phys/models/analytic_hamiltonians/fe_as_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/singleband_chain.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/analytic_hamiltonians/threeband_hubbard.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/num_traits.hpp"

namespace {

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/phys/symmetrization/";

// Double-precision tolerance for the deterministic G0 oracle
constexpr double kTolerance = 1e-12;

using Concurrency = dca::parallel::NoConcurrency;

using dca::func::dmn_0;
using dca::func::dmn_variadic;
using dca::func::function;
using BDmn = dmn_0<dca::phys::domains::electron_band_domain>;
using SDmn = dmn_0<dca::phys::domains::electron_spin_domain>;
using NuDmn = dmn_variadic<BDmn, SDmn>;

// ---------------------------------------------------------------------------------
// Case traits. THIS_TEST_DEFINES (set by CMake TEST_DEFINES) selects which one is
// instantiated; only one trait is ever live per binary, so the domain singletons are
// never shared across cases.
//
// The expectedFailing*() lists hold the ops/reps whose encoding the per-op oracle
// finds inconsistent with the deterministic G0 -- i.e. the committed red cells. They
// are filled from the harness's own first run and become the gate: any change flips the
// test. (Empty == the whole point-group is encoded correctly for this representation.)
// ---------------------------------------------------------------------------------

// CASE 1 -- square lattice with the full D4 point group (8 ops)
struct SquareD4 {
  using Scalar = double;
  using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
  static constexpr char Input[] = "square_D4_input.json";
  static constexpr int expected_num_symmetries = 8;
  // GREEN baseline: single band, no off-diagonals -- nothing for the imposition to
  // get wrong. Validates that the harness itself reports green when it should.
  static std::vector<std::string> expectedFailingReps() { return {}; }
  static std::vector<int> expectedFailingKOps() { return {}; }
  static std::vector<int> expectedFailingROps() { return {}; }
  static std::vector<int> expectedUnverifiedKOps() { return {}; }
  static constexpr int expected_num_derived_symmetries = 8;
};

// CASE 2 -- three-band Hubbard (Emery/CuO2) model on D4. This case fails the
// newly-introduced test, indicating a possible bug we expect to fix in the next
// milestone of the symmetrization project.
struct ThreebandD4 {
  using Scalar = double;
  using Lattice = dca::phys::models::ThreebandHubbard<dca::phys::domains::D4>;
  static constexpr char Input[] = "threeband_D4_input.json";
  static constexpr int expected_num_symmetries = 8;
  // The per-op cluster maps are GREEN -- each single operation's table encoding is
  // consistent with G0 -- but the production Symmetrize::execute is NOT a no-op on
  // the deterministic G0, suggesting buggy behavior in how multi-band symmetry
  // is imposed. These flip to green when the next milestone in the symmetrization
  // project rewrites the imposition.
  static std::vector<std::string> expectedFailingReps() { return {"k_iw", "r_iw", "r_tau"}; }
  static std::vector<int> expectedFailingKOps() { return {}; }
  static std::vector<int> expectedFailingROps() { return {}; }
  static std::vector<int> expectedUnverifiedKOps() { return {}; }
  static constexpr int expected_num_derived_symmetries = 8;
};

// CASE 3 -- 1D chain. Since this model has identity-only symmetry, this test exercises
// only the basic plumbing (valuable as a diagnostic).
struct SinglebandChain {
  using Scalar = double;
  // singleband_chain hardcodes no_symmetry<2> regardless of the template argument.
  using Lattice = dca::phys::models::singleband_chain<dca::phys::domains::no_symmetry<2>>;
  static constexpr char Input[] = "singleband_chain_input.json";
  static constexpr int expected_num_symmetries = 1;
  // GREEN baseline: identity-only group, single band.
  static std::vector<std::string> expectedFailingReps() { return {}; }
  static std::vector<int> expectedFailingKOps() { return {}; }
  static std::vector<int> expectedFailingROps() { return {}; }
  static std::vector<int> expectedUnverifiedKOps() { return {}; }
  static constexpr int expected_num_derived_symmetries = 4;
};

// CASE 4 -- Kagome lattice declared with no_symmetry<2>, mirroring the production
// instantiation in test/unit/phys/dca_step/cluster_solver/test_setup.hpp. This is the
// regression case for deriveOrbitalOpForOp's -1 guard: the hexagonal has the 12 D6 ops,
// but the orbital position+flavor matching in set_symmetry_matrices filters them all out
// and records (-1, -1). Deriving P from H0 recovers all 12 -- the search skips the -1
// candidate and finds the sublattice permutation geometry could not place.
struct KagomeNoSym {
  using Scalar = double;
  using Lattice = dca::phys::models::KagomeHubbard<dca::phys::domains::no_symmetry<2>>;
  static constexpr char Input[] = "kagome_input.json";
  static constexpr int expected_num_symmetries = 1;
  static std::vector<std::string> expectedFailingReps() { return {}; }
  static std::vector<int> expectedFailingKOps() { return {}; }
  static std::vector<int> expectedFailingROps() { return {}; }
  static std::vector<int> expectedUnverifiedKOps() { return {}; }
  static constexpr int expected_num_derived_symmetries = 12;
};

// CASE 5 -- FeAs (minimal two-band iron-pnictide model, PRB 77, 220503): the motivating case for
// deriving P from H0. Its two orbitals (d_xz, d_yz) share a site (a_vec = 0 for both), so
// position/flavor matching can only return the IDENTITY band image. But H0's two diagonals are the
// same dispersion with k_x and k_y exchanged, so the four D4 ops exchanging the axes (C4, C4^3, the
// two diagonal mirrors) are symmetries only with the band SWAP: geometry derives 4 ops, H0 all 8.
struct FeAs {
  using Scalar = double;
  // FeAsLattice ignores its template argument: it hardcodes FeAsPointGroup = {identity, C2}, the
  // largest group the geometry-derived permutation could support.
  using Lattice = dca::phys::models::FeAsLattice<dca::phys::domains::D4>;
  static constexpr char Input[] = "fe_as_input.json";
  static constexpr int expected_num_symmetries = 2;
  static std::vector<std::string> expectedFailingReps() { return {}; }
  static std::vector<int> expectedFailingKOps() { return {}; }
  static std::vector<int> expectedFailingROps() { return {}; }
  static std::vector<int> expectedUnverifiedKOps() { return {}; }
  static constexpr int expected_num_derived_symmetries = 8;
};

// Templated test fixture
template <typename C>
class SymmetrizeCharacterizationTest : public ::testing::Test {
protected:
  using Scalar = typename C::Scalar;
  using Lattice = typename C::Lattice;
  using Model = dca::phys::models::TightBindingModel<Lattice>;

  using Parameters =
      dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler, Model,
                                    void, dca::ClusterSolverId::CT_AUX,
                                    dca::NumericalTraits<dca::util::RealAlias<Scalar>, Scalar>>;

  using RClusterDmn = typename Parameters::RClusterDmn;
  using KClusterDmn = typename Parameters::KClusterDmn;
  using TDmn = typename Parameters::TDmn;
  using WDmn = typename Parameters::WDmn;

  using KCluster = typename KClusterDmn::parameter_type;
  using RCluster = typename RClusterDmn::parameter_type;
  using KSymDmn = typename dca::phys::domains::cluster_symmetry<KCluster>::sym_super_cell_dmn_t;
  using RSymDmn = typename dca::phys::domains::cluster_symmetry<RCluster>::sym_super_cell_dmn_t;

  SymmetrizeCharacterizationTest() : concurrency_(0, nullptr), parameters_("", concurrency_) {
    // load the per-case JSON into parameters_, derive model-level data.
    parameters_.template read_input_and_broadcast<dca::io::JSONReader>(input_dir + C::Input);
    parameters_.update_model();

    // Domain singletons are process-global and must be initialized exactly once, even
    // though gtest builds a fresh fixture per test. This static guard lets multiple
    // TYPED_TESTs share one init in a single binary -- tp_accumulator_gpu_test.cpp hits
    // the same singleton constraint but works around it by only ever running one case
    // per binary. reset() is inside the guard because only the first fixture's members
    // were constructed before update_domains() ran (so only they are mis-sized); the H0/
    // H_symmetry fills are outside because every fixture gets fresh, zero-valued members.
    static bool domain_initialized = false;
    if (!domain_initialized) {
      parameters_.update_domains();
      H0_.reset();
      H_symmetry_.reset();
      domain_initialized = true;
    }

    Parameters::model_type::initializeH0(parameters_, H0_);
    Parameters::model_type::initialize_H_symmetries(H_symmetry_);
  }

  Concurrency concurrency_;
  Parameters parameters_;

  function<std::complex<double>, dmn_variadic<NuDmn, NuDmn, KClusterDmn>> H0_;
  function<int, dmn_variadic<NuDmn, NuDmn>> H_symmetry_;
};

using TestTypes = ::testing::Types<THIS_TEST_DEFINES>;
TYPED_TEST_CASE(SymmetrizeCharacterizationTest, TestTypes);

// Largest element-wise magnitude of (a - b) over two identically-shaped functions.
template <typename F>
double maxAbsDiff(const F& a, const F& b) {
  double m = 0.;
  // std::abs correctly returns the modulus for std::complex<double>.
  for (int i = 0; i < a.size(); ++i)
    m = std::max(m, std::abs(a(i) - b(i)));
  return m;
}

// Runs fn and returns the what() string of any std::exception it throws, or "" if it ran
// to completion. The rejection fixtures use this to assert not just that the sign solver
// throws, but that it throws for the intended reason (the three branches catch three
// distinct corruptions and share two exception types, so the message is what distinguishes
// them).
template <typename F>
std::string captureThrowMessage(F&& fn) {
  try {
    fn();
  }
  catch (const std::exception& e) {
    return e.what();
  }
  return "";
}

// Location of an off-diagonal H0 coupling: the band pair (b0, b1) at momentum k, and its magnitude.
struct OffDiagTarget {
  int b0, b1, k;
  double mag;
};

// Finds the strongest off-diagonal coupling H0(b0, b1, k) (b0 != b1) over all momenta. Shared by the
// magnitude-mismatch and non-unit-ratio rejection fixtures, which perturb that one entry (and its
// Hermitian partner) in different ways.
template <typename F>
OffDiagTarget findStrongestOffDiagonal(const F& H0, int nb, int nk) {
  OffDiagTarget t{-1, -1, -1, 0.};
  for (int k = 0; k < nk; ++k)
    for (int b0 = 0; b0 < nb; ++b0)
      for (int b1 = 0; b1 < nb; ++b1)
        if (b0 != b1 && std::abs(H0(b0, 0, b1, 0, k)) > t.mag)
          t = {b0, b1, k, std::abs(H0(b0, 0, b1, 0, k))};
  return t;
}

// TEST 1: DomainActivation -- sanity gate. Reads the op counts off the symmetry
// domains and asserts they match the trait's expected_num_symmetries. Without this, a
// misconfigured case could run green simply because it had only the identity op to check.
TYPED_TEST(SymmetrizeCharacterizationTest, DomainActivation) {
  using Fixture = SymmetrizeCharacterizationTest<TypeParam>;
  const int n_k_ops = Fixture::KSymDmn::dmn_size();
  const int n_r_ops = Fixture::RSymDmn::dmn_size();

  std::cout << "[case] " << TypeParam::Input << ": cluster sym ops = " << n_k_ops
            << " (momentum), " << n_r_ops << " (real space); bands = " << BDmn::dmn_size()
            << ", k-points = " << Fixture::KClusterDmn::dmn_size() << "\n";

  EXPECT_EQ(n_k_ops, TypeParam::expected_num_symmetries);
  EXPECT_EQ(n_r_ops, TypeParam::expected_num_symmetries);
}

// TEST 2: CoarseNoOp -- the full production pipeline, asserted at 1e-12. Covers the
// four representations the existing test exercises (k/r x w/t). For each, run
// Symmetrize::execute on a copy of G0 and diff against the original. The set of
// representations that are not left invariant is compared to the expected failing set
// declared in the case trait, so all cases stay green now while the known-red cells gate
// the next milestone of the symmetrization project.
//
// Unlike the per-op maps, this runs the whole production pipeline (spin + cluster +
// time/frequency channels and the averaging/norm), so it catches defects in how the
// operations are combined.
TYPED_TEST(SymmetrizeCharacterizationTest, CoarseNoOp) {
  using Fixture = SymmetrizeCharacterizationTest<TypeParam>;
  using Parameters = typename Fixture::Parameters;
  using KClusterDmn = typename Fixture::KClusterDmn;
  using RClusterDmn = typename Fixture::RClusterDmn;
  using TDmn = typename Fixture::TDmn;
  using WDmn = typename Fixture::WDmn;

  // Symmetrize a copy with the production code; return max element-wise change.
  auto noOpResidual = [](auto G0) {
    auto sym = G0;
    dca::phys::Symmetrize<Parameters>::execute(sym, false);
    return maxAbsDiff(G0, sym);
  };

  function<std::complex<double>, dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>> G0_k_w;
  dca::phys::compute_G0_k_w(this->H0_, this->parameters_.get_chemical_potential(), 1, G0_k_w);

  function<std::complex<double>, dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>> G0_r_w;
  dca::math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(G0_k_w, G0_r_w);

  function<std::complex<double>, dmn_variadic<NuDmn, NuDmn, KClusterDmn, TDmn>> G0_k_t;
  dca::phys::compute_G0_k_t(this->H0_, this->parameters_.get_chemical_potential(),
                            this->parameters_.get_beta(), G0_k_t);

  function<std::complex<double>, dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>> G0_r_t_cmplx;
  dca::math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(G0_k_t, G0_r_t_cmplx);
  function<double, dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>> G0_r_t;
  G0_r_t = dca::func::util::real(G0_r_t_cmplx, dca::ImagCheck::FAIL);

  // Run the oracle on all four reps and collect (name, residual) pairs
  const std::vector<std::pair<std::string, double>> residuals = {
      {"k_iw", noOpResidual(G0_k_w)},
      {"r_iw", noOpResidual(G0_r_w)},
      {"k_tau", noOpResidual(G0_k_t)},
      {"r_tau", noOpResidual(G0_r_t)}};

  // Turn residuals into a red/green list and diff against the expected failing set
  std::vector<std::string> failing_reps;
  std::cout << "[coarse no-op residuals] " << TypeParam::Input << "\n";
  for (const auto& [rep, residual] : residuals) {
    const bool failed = residual >= kTolerance;
    std::cout << "  " << rep << ": " << (failed ? "FAIL" : "PASS")
              << "  max|G0 - sym(G0)| = " << residual << "\n";
    if (failed)
      failing_reps.push_back(rep);
  }

  // Compare as a set: sort both sides so the assertion does not depend on the residuals
  // iteration order (which carries no meaning).
  std::sort(failing_reps.begin(), failing_reps.end());
  auto expected_reps = TypeParam::expectedFailingReps();
  std::sort(expected_reps.begin(), expected_reps.end());
  EXPECT_EQ(failing_reps, expected_reps)
      << "Coarse no-op red/green map changed vs the expected failing set.";
}

// TEST 3: PerOpMapMomentum -- per-operation map (momentum channel). For each
// operation s, mirror the momentum band-aware executeCluster relation on G0(k, iw):
//   G0(b0, b1, k) == sign_s * G0(b0_new, b1_new, k_new)   when sign_s != 0,
// where (k_new, b*_new) come from get_symmetry_matrix() and sign_s from
// transformationSignOfK. The set of ops with a real violation is compared to
// the expected failing set.
TYPED_TEST(SymmetrizeCharacterizationTest, PerOpMapMomentum) {
  using Fixture = SymmetrizeCharacterizationTest<TypeParam>;
  using KClusterDmn = typename Fixture::KClusterDmn;
  using WDmn = typename Fixture::WDmn;
  using KCluster = typename Fixture::KCluster;
  using Lattice = typename Fixture::Lattice;

  // Build the (k, iw) oracle exactly as in CoarseNoOp.
  function<std::complex<double>, dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>> G0;
  dca::phys::compute_G0_k_w(this->H0_, this->parameters_.get_chemical_potential(), 1, G0);

  const auto& sym = dca::phys::domains::cluster_symmetry<KCluster>::get_symmetry_matrix();
  const int n_ops = Fixture::KSymDmn::dmn_size();
  const int nb = BDmn::dmn_size();
  const int ns = SDmn::dmn_size();
  const int nk = KClusterDmn::dmn_size();
  const int nw = WDmn::dmn_size();

  std::vector<int> failing_ops;
  std::cout << "[per-op momentum map] " << TypeParam::Input << "\n";
  // ---- Outer loop over symmetry ops s. Each iteration is one independent check.
  for (int s = 0; s < n_ops; ++s) {
    double max_viol = 0.;
    bool any_checked = false;  // Tracks whether op s constrained anything (else N/A)
    // ---- Sweep every (k, b0, b1) band-pair at this op
    for (int k = 0; k < nk; ++k)
      for (int b0 = 0; b0 < nb; ++b0)
        for (int b1 = 0; b1 < nb; ++b1) {
          const double sign = Lattice::transformationSignOfK(b0, b1, s);
          if (sign == 0)  // op does not constrain this band pair -> N/A
            continue;
          const int k_new = sym(k, b0, s).first;  // mirrors production ("FIXME: b0 -> b1")
          const int b0_new = sym(k, b0, s).second;
          const int b1_new = sym(k, b1, s).second;
          // ---- Invariance check, swept over spin and frequency. Asserts
          // G0(b0,b1,k,w) == sign * G0(b0_new,b1_new,k_new,w). Accumulate worst
          // violation rather than failing on first, so the printed number is the
          // real max error for this op.
          for (int sp = 0; sp < ns; ++sp)
            for (int w = 0; w < nw; ++w) {
              any_checked = true;
              const auto lhs = G0(b0, sp, b1, sp, k, w);
              const auto rhs = sign * G0(b0_new, sp, b1_new, sp, k_new, w);
              max_viol = std::max(max_viol, std::abs(lhs - rhs));
            }
        }
    // Classify this op: FAIL only if it actually checked something AND exceeded tolerance
    const bool failed = any_checked && max_viol >= kTolerance;
    std::cout << "  op " << s << ": " << (any_checked ? (failed ? "FAIL" : "PASS") : "N/A ")
              << "  max|G0 - S.G0| = " << max_viol << "\n";
    if (failed)
      failing_ops.push_back(s);
  }

  std::sort(failing_ops.begin(), failing_ops.end());
  auto expected_ops = TypeParam::expectedFailingKOps();
  std::sort(expected_ops.begin(), expected_ops.end());
  EXPECT_EQ(failing_ops, expected_ops)
      << "Momentum per-op red/green map changed vs the expected failing set.";
}

// TEST 4: PerOpMapMomentumUS -- the shadow orbital-op record, derived from H0. Solve each U_S as
// a signed permutation directly from H0(k) by enumerating the sign assignments consistent with the
// fold-corrected H0 couplings (no transformationSignOf*), verify it is a well-formed unitary signed
// permutation, and re-run the per-op momentum no-op as the matrix conjugation
// G0(k) == (U_S V) G0(k_new) (U_S V)^dagger, where V = diag(fold_phase(k, ., s)) redresses the
// intrinsic U_S for the folded representatives the table stores -- rather than the element-wise
// scalar relation of TEST 3. It must reproduce TEST 3's red/green map; its purpose is to prove the
// consolidated, H0-derived matrix record is faithful and correct, not to find new bugs.
TYPED_TEST(SymmetrizeCharacterizationTest, PerOpMapMomentumUS) {
  using Fixture = SymmetrizeCharacterizationTest<TypeParam>;
  using KClusterDmn = typename Fixture::KClusterDmn;
  using WDmn = typename Fixture::WDmn;
  using KCluster = typename Fixture::KCluster;

  // Populate the shadow orbital-op table from H0 alone; throws if the op is not a signed-permutation
  // symmetry of H0.
  dca::phys::solveOrbitalOpSignsFromH0<KCluster>(this->H0_);
  const auto& u_s = dca::phys::domains::cluster_symmetry<KCluster>::get_orbital_op();
  const auto& sym = dca::phys::domains::cluster_symmetry<KCluster>::get_symmetry_matrix();
  const auto& fold_phase = dca::phys::domains::cluster_symmetry<KCluster>::get_fold_phase();

  const int n_ops = Fixture::KSymDmn::dmn_size();
  const int nb = BDmn::dmn_size();
  const int ns = SDmn::dmn_size();
  const int nk = KClusterDmn::dmn_size();
  const int nw = WDmn::dmn_size();

  // ---- Check 1 (structural): each U_S is a signed permutation -- exactly one +/-1 per row and column.
  // (A real matrix with one unit-magnitude entry per row and per column is orthogonal by
  // construction, so there is no separate U_S U_S^T = I check to run.)
  for (int s = 0; s < n_ops; ++s)
    for (int r = 0; r < nb; ++r) {
      int row_nnz = 0, col_nnz = 0;
      for (int c = 0; c < nb; ++c) {
        if (u_s(r, c, s) != 0.) {
          ++row_nnz;
          EXPECT_DOUBLE_EQ(std::abs(u_s(r, c, s)), 1.0);
        }
        if (u_s(c, r, s) != 0.)
          ++col_nnz;
      }
      EXPECT_EQ(row_nnz, 1) << "U_S op " << s << " row " << r << " is not a permutation row";
      EXPECT_EQ(col_nnz, 1) << "U_S op " << s << " col " << r << " is not a permutation column";
    }

  // ---- Check 2 (action): per-op no-op via conjugation  G0(k) == (U_S V) G0(k_new) (U_S V)^dagger
  // (all real). U_S is the intrinsic orbital transform; V = diag(fold_phase(k, ., s)) restores the
  // folding phase of the representative k_new = S k - G the table stores, at the image (contracted)
  // band indices.
  function<std::complex<double>, dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>> G0;
  dca::phys::compute_G0_k_w(this->H0_, this->parameters_.get_chemical_potential(), 1, G0);

  std::vector<int> failing_ops;
  std::cout << "[per-op momentum map, U_S conjugation] " << TypeParam::Input << "\n";
  for (int s = 0; s < n_ops; ++s) {
    double max_viol = 0.;
    for (int k = 0; k < nk; ++k) {
      const int k_new = sym(k, 0, s).first;  // k-image is band-independent
      for (int sp = 0; sp < ns; ++sp)
        for (int w = 0; w < nw; ++w)
          for (int b0 = 0; b0 < nb; ++b0)
            for (int b1 = 0; b1 < nb; ++b1) {
              std::complex<double> conj_val = 0.;
              for (int a = 0; a < nb; ++a)
                for (int b = 0; b < nb; ++b)
                  conj_val += u_s(b0, a, s) * fold_phase(k, a, s) * G0(a, sp, b, sp, k_new, w) *
                              fold_phase(k, b, s) * u_s(b1, b, s);
              max_viol = std::max(max_viol, std::abs(G0(b0, sp, b1, sp, k, w) - conj_val));
            }
    }
    const bool failed = max_viol >= kTolerance;
    std::cout << "  op " << s << ": " << (failed ? "FAIL" : "PASS")
              << "  max|G0 - U_S.G0.U_S^dag| = " << max_viol << "\n";
    if (failed)
      failing_ops.push_back(s);
  }

  std::sort(failing_ops.begin(), failing_ops.end());
  auto expected_ops = TypeParam::expectedFailingKOps();
  std::sort(expected_ops.begin(), expected_ops.end());
  EXPECT_EQ(failing_ops, expected_ops)
      << "U_S conjugation per-op momentum map changed vs the expected failing set.";
}

// TEST 5: GroupShadowCrossCheck -- the headline U_S verification statement. The candidate ops
// in the symmetry table come from the declared point group (DCA_point_group), filtered by lattice
// geometry; DomainActivation already pins their count to the declared order. This test closes the
// loop at the Hamiltonian level: the set of ops the sign-consistency solve accepts (the H0-verified
// group) must equal the full declared/realized group. Equivalently, no declared op may be merely a
// lattice symmetry that H0 does not respect. The set of rejected ops is compared to the per-case
// expectation (empty on every model included in this test). When candidate ops are derived from H0,
// this same predicate becomes the actual filter that carves the real group out of a candidate table;
// here it is validated against the declared group as a known-correct oracle.
TYPED_TEST(SymmetrizeCharacterizationTest, GroupShadowCrossCheck) {
  using Fixture = SymmetrizeCharacterizationTest<TypeParam>;
  using KCluster = typename Fixture::KCluster;
  const int n_ops = Fixture::KSymDmn::dmn_size();

  // The H0-verified group: ops the sign-consistency solve accepts (non-throwing).
  const std::vector<int> verified = dca::phys::verifiedSymmetryOps<KCluster>(this->H0_);

  // Anything the table declares but H0 does not verify.
  std::vector<int> rejected;
  for (int s = 0; s < n_ops; ++s)
    if (std::find(verified.begin(), verified.end(), s) == verified.end())
      rejected.push_back(s);

  std::cout << "[group shadow] " << TypeParam::Input << ": " << verified.size() << "/" << n_ops
            << " declared ops H0-verified\n";

  std::sort(rejected.begin(), rejected.end());
  auto expected_rejected = TypeParam::expectedUnverifiedKOps();
  std::sort(expected_rejected.begin(), expected_rejected.end());
  EXPECT_EQ(rejected, expected_rejected)
      << "H0-verified group differs from the declared group: the set of declared ops the "
         "sign-consistency check rejects changed vs the expectation.";
}

// TEST 5c: DerivedGroupReport -- the production derive-and-report step. update_domains now derives
// each 2D model's point group from scratch (2D holohedry pool -> geometry filter -> H0 gate),
// reports any divergence from the declared group, and keeps symmetrizing with the declared group.
// This exercises the same routine directly and asserts the three things production relies on:
// (1) the declared symmetry state is restored, so the derivation is invisible downstream
// (2) the report is silent when derived == declared, and names the divergence otherwise
// (3) the derived-group size matches the trait expectation
TYPED_TEST(SymmetrizeCharacterizationTest, DerivedGroupReport) {
  using Fixture = SymmetrizeCharacterizationTest<TypeParam>;
  using RClusterDmn = typename Fixture::RClusterDmn;
  using KCluster = typename Fixture::KCluster;
  using Model = typename Fixture::Model;
  using CS_k = dca::phys::domains::cluster_symmetry<KCluster>;

  const int n_ops_before = Fixture::KSymDmn::dmn_size();
  const auto symmetry_matrix_before = CS_k::get_symmetry_matrix();

  const std::string report =
      dca::phys::deriveAndComparePointGroup<RClusterDmn, Model>(this->parameters_);

  // (1) restoration: production must be left in exactly the declared state.
  ASSERT_EQ(Fixture::KSymDmn::dmn_size(), n_ops_before);
  const auto& symmetry_matrix_after = CS_k::get_symmetry_matrix();
  for (int i = 0; i < symmetry_matrix_before.size(); ++i)
    ASSERT_EQ(symmetry_matrix_after(i), symmetry_matrix_before(i))
        << "symmetry table changed at linear index " << i;

  // (2) + (3): silent on match; on divergence the report names the derived-group size and the
  // under-declaration.
  if (TypeParam::expected_num_derived_symmetries == TypeParam::expected_num_symmetries) {
    EXPECT_TRUE(report.empty()) << "unexpected divergence report:\n" << report;
  }
  else {
    const std::string expected_size = "derived group: " +
                                      std::to_string(TypeParam::expected_num_derived_symmetries) +
                                      " op(s)";
    EXPECT_NE(report.find(expected_size), std::string::npos) << report;
    EXPECT_NE(report.find("under-declared"), std::string::npos) << report;
    std::cout << "[derived group] " << TypeParam::Input << ":\n" << report;
  }
}

// TEST 6: MappedPointShadow -- the mapped_point accessor. The band-independent ".first" of the
// legacy symmetry matrix (the image of each momentum under each op) is promoted into its own
// (cluster-point x op) accessor. Verify it reproduces .first for every band, which confirms both
// that the shadow accessor is faithful and that the momentum image really is band-independent.
TYPED_TEST(SymmetrizeCharacterizationTest, MappedPointShadow) {
  using Fixture = SymmetrizeCharacterizationTest<TypeParam>;
  using KCluster = typename Fixture::KCluster;
  using KClusterDmn = typename Fixture::KClusterDmn;

  const auto& sym = dca::phys::domains::cluster_symmetry<KCluster>::get_symmetry_matrix();
  const auto& mapped_point = dca::phys::domains::cluster_symmetry<KCluster>::get_mapped_point();

  const int n_ops = Fixture::KSymDmn::dmn_size();
  const int nb = BDmn::dmn_size();
  const int nk = KClusterDmn::dmn_size();

  for (int s = 0; s < n_ops; ++s)
    for (int k = 0; k < nk; ++k)
      for (int b = 0; b < nb; ++b)
        EXPECT_EQ(mapped_point(k, s), sym(k, b, s).first)
            << "mapped_point disagrees with legacy .first at (k=" << k << ", b=" << b << ", op=" << s
            << ")";
}

// TEST 7: PerOpMapRealSpace -- mirror of TEST 3 but in real space / tau, and with one
// key difference: it skips off-diagonal band pairs (b0!=b1) because production's real-space
// executeCluster does the same. That skip is itself one of the suspected bugs, so this test
// characterizes production's behavior, off-diag blind spot included.
TYPED_TEST(SymmetrizeCharacterizationTest, PerOpMapRealSpace) {
  using Fixture = SymmetrizeCharacterizationTest<TypeParam>;
  using KClusterDmn = typename Fixture::KClusterDmn;
  using RClusterDmn = typename Fixture::RClusterDmn;
  using TDmn = typename Fixture::TDmn;
  using RCluster = typename Fixture::RCluster;
  using Lattice = typename Fixture::Lattice;

  // Build (k,tau), FT to (r,tau), drop the imaginary part -> real G0(r,tau).
  function<std::complex<double>, dmn_variadic<NuDmn, NuDmn, KClusterDmn, TDmn>> G0_k_t;
  dca::phys::compute_G0_k_t(this->H0_, this->parameters_.get_chemical_potential(),
                            this->parameters_.get_beta(), G0_k_t);
  function<std::complex<double>, dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>> G0_cmplx;
  dca::math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(G0_k_t, G0_cmplx);
  function<double, dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>> G0;
  G0 = dca::func::util::real(G0_cmplx, dca::ImagCheck::FAIL);

  const auto& sym = dca::phys::domains::cluster_symmetry<RCluster>::get_symmetry_matrix();
  const int n_ops = Fixture::RSymDmn::dmn_size();
  const int nb = BDmn::dmn_size();
  const int ns = SDmn::dmn_size();
  const int nr = RClusterDmn::dmn_size();
  const int nt = TDmn::dmn_size();

  std::vector<int> failing_ops;
  std::cout << "[per-op real-space map] " << TypeParam::Input << "\n";
  for (int s = 0; s < n_ops; ++s) {
    double max_viol = 0.;
    bool any_checked = false;
    for (int r = 0; r < nr; ++r)
      for (int b0 = 0; b0 < nb; ++b0)
        for (int b1 = 0; b1 < nb; ++b1) {
          if (b0 != b1)  // off-diagonal skip in executeCluster -> N/A
            continue;
          // Real-space phase for this diagonal element under op s; 0 => N/A.
          const double sign = Lattice::transformationSignOfR(b0, b1, s);
          if (sign == 0)
            continue;
          const int r_new = sym(r, 0, s).first;
          const int b0_new = sym(r, b0, s).second;
          const int b1_new = sym(0, b1, s).second;
          // Invariance check over spin + imaginary time, real arithmetic this time.
          for (int sp = 0; sp < ns; ++sp)
            for (int t = 0; t < nt; ++t) {
              any_checked = true;
              const double lhs = G0(b0, sp, b1, sp, r, t);
              const double rhs = sign * G0(b0_new, sp, b1_new, sp, r_new, t);
              max_viol = std::max(max_viol, std::abs(lhs - rhs));
            }
        }
    const bool failed = any_checked && max_viol >= kTolerance;
    std::cout << "  op " << s << ": " << (any_checked ? (failed ? "FAIL" : "PASS") : "N/A ")
              << "  max|G0 - S.G0| = " << max_viol << "\n";
    if (failed)
      failing_ops.push_back(s);
  }

  std::sort(failing_ops.begin(), failing_ops.end());
  auto expected_ops = TypeParam::expectedFailingROps();
  std::sort(expected_ops.begin(), expected_ops.end());
  EXPECT_EQ(failing_ops, expected_ops)
      << "Real-space per-op red/green map changed vs the expected failing set.";
}

// TEST 8 (a/b/c/d): rejection fixtures for the sign-consistency solver.
//
// solveOrbitalOpSignsFromH0 derives the permutation from H0, throwing only when no band
// permutation reproduces it. No shipped model triggers that on a correct H0 -- the spine cases all
// pass -- so we perturb a copy of the (exactly symmetric) H0 (8a-8c) to break one op past every
// permutation, or corrupt the table (8d) to show the search recovers what geometry cannot place.
//
// These need a multi-orbital model. A single band has no off-diagonal coupling to corrupt,
// and on the spine's square cluster every nonzero single-band entry sits at a k that is a
// fixed point of every op (the BZ center and corner), so no single-entry perturbation can
// break a symmetry there. Threeband (d, px, py) is the driver; square and singleband skip.

// 8a/8b -- two ways a single broken off-diagonal coupling breaks the op past any permutation. Both
// perturb the strongest off-diagonal H0 entry (and its Hermitian partner), differing only in how:
//   * vanish it (a coupling present on one side of the relation but gone on the other);
//   * scale it by 1.5 (the H0 ratio that would fix a sign product is no longer a real +/-1).
// Either way the search exhausts every permutation and throws the unified message.
TYPED_TEST(SymmetrizeCharacterizationTest, RejectsBrokenOffDiagonalCoupling) {
  using Fixture = SymmetrizeCharacterizationTest<TypeParam>;
  using KCluster = typename Fixture::KCluster;
  const int nb = BDmn::dmn_size();
  const int nk = Fixture::KClusterDmn::dmn_size();
  if (nb < 2) {
    std::cout << "[skip] " << TypeParam::Input
              << ": needs a multi-orbital model with off-diagonal couplings.\n";
    return;
  }

  const auto t = findStrongestOffDiagonal(this->H0_, nb, nk);
  ASSERT_GT(t.mag, kTolerance) << "no off-diagonal coupling to perturb.";

  // The perturbation is only visible to an op that maps the corrupted entry to a *different*
  // (unperturbed) H0 entry. Both the entry and its Hermitian partner are perturbed, so an op
  // landing on either sees a consistent H0 and cannot throw. An identity-only group (e.g. a
  // multi-band lattice declared with no_symmetry) therefore has nothing to detect: skip.
  {
    const auto& sym = dca::phys::domains::cluster_symmetry<KCluster>::get_symmetry_matrix();
    const int n_ops = Fixture::KSymDmn::dmn_size();
    bool movable = false;
    for (int s = 0; s < n_ops && !movable; ++s) {
      const int k_new = sym(t.k, 0, s).first;
      const int b0_new = sym(0, t.b0, s).second;
      const int b1_new = sym(0, t.b1, s).second;
      const bool fixes = k_new == t.k && b0_new == t.b0 && b1_new == t.b1;
      const bool hermitian_partner = k_new == t.k && b0_new == t.b1 && b1_new == t.b0;
      movable = !fixes && !hermitian_partner;
    }
    if (!movable) {
      std::cout << "[skip] " << TypeParam::Input
                << ": no op moves the corrupted coupling to an unperturbed entry.\n";
      return;
    }
  }

  // Vanish: the op connecting this entry to its (still nonzero) image sees it gone on one side only.
  {
    auto H0p = this->H0_;
    H0p(t.b0, 0, t.b1, 0, t.k) = 0.;
    H0p(t.b1, 0, t.b0, 0, t.k) = 0.;
    const std::string msg =
        captureThrowMessage([&] { dca::phys::solveOrbitalOpSignsFromH0<KCluster>(H0p); });
    ASSERT_FALSE(msg.empty()) << "sign solver did not reject a vanished coupling.";
    EXPECT_NE(msg.find("for any band permutation"), std::string::npos) << "actual: " << msg;
  }

  // Scale by 1.5: the entry stays nonzero but its ratio to the image is not +/-1.
  {
    auto H0p = this->H0_;
    H0p(t.b0, 0, t.b1, 0, t.k) *= 1.5;
    H0p(t.b1, 0, t.b0, 0, t.k) *= 1.5;
    const std::string msg =
        captureThrowMessage([&] { dca::phys::solveOrbitalOpSignsFromH0<KCluster>(H0p); });
    ASSERT_FALSE(msg.empty()) << "sign solver did not reject a non-unit coupling ratio.";
    EXPECT_NE(msg.find("for any band permutation"), std::string::npos) << "actual: " << msg;
  }
}

// 8c -- inconsistent signs: pairwise sign products that are each a valid +/-1 but cannot be
// jointly satisfied by any assignment. At the zone corner all three bands couple (d-px, d-py,
// px-py), so the d-px and d-py signs force the px-py sign product, and a px<->py swap op
// reads the px-py sign off H0 directly. Flipping one direction of the px-py coupling
// makes the direct reading contradict the forced product: no sign vector works -- and, because the
// flip is asymmetric, no other band permutation rescues it either, so the search still throws.
TYPED_TEST(SymmetrizeCharacterizationTest, RejectsInconsistentSigns) {
  using Fixture = SymmetrizeCharacterizationTest<TypeParam>;
  using KCluster = typename Fixture::KCluster;
  const int nb = BDmn::dmn_size();
  const int nk = Fixture::KClusterDmn::dmn_size();
  if (nb < 3) {
    std::cout << "[skip] " << TypeParam::Input << ": needs >=3 coupled bands to form a sign cycle.\n";
    return;
  }

  auto H0p = this->H0_;
  // The zone corner is the only momentum where the px-py coupling H0(1,2) is nonzero (both
  // sin(k/2) factors nonzero); there all three bands couple.
  int kc = -1;
  for (int k = 0; k < nk; ++k)
    if (std::abs(H0p(1, 0, 2, 0, k)) > kTolerance) {
      kc = k;
      break;
    }
  ASSERT_GE(kc, 0) << "no zone-corner px-py coupling found.";

  // An op that swaps bands 1<->2 and fixes the corner reads the px-py sign as the raw H0
  // ratio; without one the contradiction never arises.
  const auto& sym = dca::phys::domains::cluster_symmetry<KCluster>::get_symmetry_matrix();
  const int n_ops = Fixture::KSymDmn::dmn_size();
  bool have_swap = false;
  for (int s = 0; s < n_ops; ++s)
    if (sym(0, 1, s).second == 2 && sym(0, 2, s).second == 1 && sym(kc, 0, s).first == kc) {
      have_swap = true;
      break;
    }
  if (!have_swap) {
    std::cout << "[skip] " << TypeParam::Input << ": no px<->py swap op fixing the zone corner.\n";
    return;
  }

  H0p(1, 0, 2, 0, kc) *= -1.;  // one direction only, on purpose.

  const std::string msg =
      captureThrowMessage([&] { dca::phys::solveOrbitalOpSignsFromH0<KCluster>(H0p); });
  ASSERT_FALSE(msg.empty()) << "sign solver did not reject an inconsistent sign system.";
  EXPECT_NE(msg.find("for any band permutation"), std::string::npos) << "actual: " << msg;
}

// TEST 9: the folding phase must not be absorbed into U_S (PR #368 review). At a zone-boundary
// momentum the folded representative of the mirror image coincides with the original point, so the
// raw H0 ratio is +1 there; only the fold correction phi = e^{i G.a_b} exposes the intrinsic odd
// parity of p_x. The scenario is the threeband d-p_x block at k = (pi, 0) under the k_x -> -k_x
// mirror: the image (-pi, 0) folds back to (pi, 0) with G = (2 pi, 0), and a_px = (1/2, 0) gives
// phi_px = cos(pi) = -1 while phi_d = +1. Without the correction the solve returns sigma_px = +1
// -- the representative-convention-dressed sign -- instead of the intrinsic -1. Mock-based and
// case-independent (unlike the fixtures above): the dressing is uniform across every momentum of
// the shipped clusters, so no perturbation of a real H0 can expose it; the single boundary point
// with an explicit fold functor is the minimal reproducer.
TEST(SolveOrbitalOpSigns, BoundaryFoldDoesNotMaskPxParity) {
  // One momentum, one op: the mirror image folds back onto the same representative, no band
  // permutation.
  struct Sym {
    std::pair<int, int> operator()(int /*k*/, int b, int /*s*/) const {
      return {0, b};
    }
  };
  // The d-p_x block of the threeband H0 at k = (pi, 0) with t_pd = 1: 2 i t_pd sin(k_x / 2).
  struct H0 {
    std::complex<double> operator()(int b0, int /*s0*/, int b1, int /*s1*/, int /*k*/) const {
      const std::complex<double> i(0., 1.);
      if (b0 == 0 && b1 == 1)
        return 2. * i;
      if (b0 == 1 && b1 == 0)
        return -2. * i;  // Hermitian partner.
      return 0.;
    }
  };
  // phi_b = cos(G . a_b) with G = (2 pi, 0): -1 for p_x (band 1, a = (1/2, 0)), +1 for d (at the
  // origin) and p_y (a = (0, 1/2), orthogonal to G).
  struct Fold {
    double operator()(int /*k*/, int b, int /*s*/) const {
      return b == 1 ? -1. : 1.;
    }
  };
  // Output stand-in for cluster_symmetry::get_orbital_op(): the solver writes the solved U_S
  // here (band row, band column; the single op collapses the op index).
  struct USTable {
    double data[3][3] = {};
    double& operator()(int r, int c, int /*s*/) {
      return data[r][c];
    }
  };

  USTable u_s;
  dca::phys::detail::deriveOrbitalOpForOp(0, 3, 1, Sym{}, Fold{}, H0{}, u_s);

  EXPECT_DOUBLE_EQ(u_s(0, 0, 0), 1.);
  EXPECT_DOUBLE_EQ(u_s(1, 1, 0), -1.);  // The intrinsic p_x parity, not the fold-dressed +1.
  EXPECT_DOUBLE_EQ(u_s(2, 2, 0), 1.);
}

// 8d -- recovery from a missing band image: set_symmetry_matrices records -1 when its position/
// flavor matching finds no admissible image. Unlike 8a-8c this perturbs the symmetry table, not H0:
// the solver skips the -1 candidate (no OOB) and recovers the op from the search, so the populator
// succeeds and the classifier still includes it. The table is static state, so it gets restored.
TYPED_TEST(SymmetrizeCharacterizationTest, RecoversFromMissingBandImage) {
  using Fixture = SymmetrizeCharacterizationTest<TypeParam>;
  using KCluster = typename Fixture::KCluster;

  const std::vector<int> baseline = dca::phys::verifiedSymmetryOps<KCluster>(this->H0_);
  ASSERT_FALSE(baseline.empty()) << "no verified op to corrupt (identity should always pass).";
  const int s_target = baseline.front();

  auto& sym = dca::phys::domains::cluster_symmetry<KCluster>::get_symmetry_matrix();
  const int saved_image = sym(0, 0, s_target).second;
  sym(0, 0, s_target).second = -1;

  // The populator must NOT throw: the -1 makes the geometric candidate inadmissible, but the search
  // should continue
  const std::string msg =
      captureThrowMessage([&] { dca::phys::solveOrbitalOpSignsFromH0<KCluster>(this->H0_); });
  const std::vector<int> verified = dca::phys::verifiedSymmetryOps<KCluster>(this->H0_);

  sym(0, 0, s_target).second = saved_image;

  EXPECT_TRUE(msg.empty()) << "solver failed to recover from a -1 valued band image: " << msg;
  EXPECT_EQ(verified, baseline)
      << "the H0 search did not recover the op whose geometric band image was corrupted.";
}

}  // namespace
