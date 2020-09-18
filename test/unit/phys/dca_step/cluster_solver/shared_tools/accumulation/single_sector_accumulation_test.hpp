// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file provides a common setup for accumulation tests. It computes a mock configuration and
// M matrix over a single spin sector.

#ifndef TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_SINGLE_SECTOR_ACCUMULATION_TEST_HPP
#define TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_SINGLE_SECTOR_ACCUMULATION_TEST_HPP

#include <cmath>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"

namespace dca {
namespace testing {

class PositiveFrq {
public:
  using element_type = double;
  static std::size_t get_size() {
    return size_;
  }
  static void initialize(const int n_frq) {
    size_ = n_frq;
  }

private:
  static int size_;
};
int PositiveFrq::size_ = -1;

template <int size>
class MockClusterDmn {
public:
  using element_type = int;

  static int get_size() {
    return size;
  }
  static int subtract(int i, int j) {
    if (i >= size || j >= size)
      throw(std::out_of_range(__FUNCTION__));
    return (j - i + size) % size;
  }

  static const auto& get_subtract_matrix() {
    auto initialize_subtract_matrix = []() {
      linalg::Matrix<int, linalg::CPU> m(size);
      for (int j = 0; j < size; ++j)
        for (int i = 0; i < size; ++i)
          m(i, j) = subtract(i, j);
      return m;
    };
    static auto sub_matrix = initialize_subtract_matrix();
    return sub_matrix;
  }

  static int add(int i, int j) {
    if (i >= size || j >= size)
      throw(std::out_of_range(__FUNCTION__));
    return (j + i) % size;
  }

  static const auto& get_add_matrix() {
    auto initialize_add_matrix = []() {
      linalg::Matrix<int, linalg::CPU> m(size);
      for (int j = 0; j < size; ++j)
        for (int i = 0; i < size; ++i)
          m(i, j) = add(i, j);
      return m;
    };

    static const auto m_add = initialize_add_matrix();
    return m_add;
  }
};

struct Vertex {
  double get_tau() const {
    return tau_;
  }
  int get_band() const {
    return b_;
  }
  int get_left_band() const {
    return b_;
  }
  int get_right_band() const {
    return b_;
  }
  int get_r_site() const {
    return r_;
  }
  int get_left_site() const {
    return r_;
  }
  int get_right_site() const {
    return r_;
  }

  int b_;
  int r_;
  double tau_;
};

namespace {
// Flag for single initialization when multiple types are used.
bool single_sector_accumulator_test_initialized = false;
}  // namespace

template <typename Real = double, int n_bands = 2, int n_sites = 3, int n_frqs = 64>
class SingleSectorAccumulationTest : public ::testing::Test {
public:
  using Complex = std::complex<Real>;

  using RDmn = dca::func::dmn_0<MockClusterDmn<n_sites>>;
  using FreqDmn = dca::func::dmn_0<dca::phys::domains::frequency_domain>;
  using PosFreqDmn = dca::func::dmn_0<PositiveFrq>;
  using BDmn = dca::func::dmn_0<dca::phys::domains::electron_band_domain>;

  using Configuration = std::vector<Vertex>;
  using Matrix = dca::linalg::Matrix<Real, dca::linalg::CPU>;

  using F_w_w =
      dca::func::function<Complex, dca::func::dmn_variadic<BDmn, BDmn, RDmn, RDmn, FreqDmn, FreqDmn>>;

  static double get_beta() {
    return beta_;
  }

public:
  struct MockParameters {
    constexpr static int bands = n_bands;
    using lattice_type = dca::phys::models::square_lattice<dca::phys::domains::no_symmetry<2>>;
  };

  static void SetUpTestCase() {
    if (!single_sector_accumulator_test_initialized) {
      // Initialize the frequency domains.
      dca::phys::domains::frequency_domain::initialize(beta_, n_frqs);
      PositiveFrq::initialize(n_frqs);
      // Initialize the band domain.
      BDmn::parameter_type::initialize(MockParameters());

      single_sector_accumulator_test_initialized = true;
    }
  }

  void SetUp() {}

  auto compute2DFTBaseline() const -> F_w_w;

  void prepareConfiguration(Configuration& config, Matrix& M, const int n) const {
    prepareConfiguration(config, M, BDmn::dmn_size(), RDmn::dmn_size(), get_beta(), n);
  }

  static void prepareConfiguration(Configuration& config, Matrix& M, const int nb, const int nr,
                                   const double beta, const int n);

protected:
  static constexpr double beta_ = 2.;
  Configuration configuration_;
  Matrix M_;
};

template <typename Real, int n_bands, int n_sites, int n_frqs>
void SingleSectorAccumulationTest<Real, n_bands, n_sites, n_frqs>::prepareConfiguration(
    Configuration& config, Matrix& M, const int nb, const int nr, const double beta, const int n) {
  config.resize(n);
  M.resize(n);
  static dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 0);

  for (int i = 0; i < n; ++i) {
    const double tau = rng() * beta;
    const int r = rng() * nr;
    const int b = rng() * nb;
    config[i] = Vertex{b, r, tau};
  }

  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) {
      M(i, j) = 2 * rng() - 1.;
    }
}

template <typename Real, int n_bands, int n_sites, int n_frqs>
auto SingleSectorAccumulationTest<Real, n_bands, n_sites, n_frqs>::compute2DFTBaseline() const
    -> F_w_w {
  F_w_w f_w("2D frequency transform baseline.");
  const std::complex<Real> imag(0, 1);

  for (int w_ind2 = 0; w_ind2 < FreqDmn::dmn_size(); ++w_ind2) {
    const Real w_val2 = FreqDmn::get_elements()[w_ind2];
    for (int w_ind1 = 0; w_ind1 < FreqDmn::dmn_size(); ++w_ind1) {
      const Real w_val1 = FreqDmn::get_elements()[w_ind1];
      for (int j = 0; j < configuration_.size(); ++j) {
        const Real t_val2 = configuration_[j].get_tau();
        const int b2 = configuration_[j].b_;
        const int r2 = configuration_[j].r_;
        for (int i = 0; i < configuration_.size(); ++i) {
          const Real t_val1 = configuration_[i].get_tau();
          const int b1 = configuration_[i].b_;
          const int r1 = configuration_[i].r_;

          const auto f_t = M_(i, j);
          f_w(b1, b2, r1, r2, w_ind1, w_ind2) +=
              f_t * std::exp(imag * (t_val1 * w_val1 - t_val2 * w_val2));
        }
      }
    }
  }

  return f_w;
}

}  // namespace testing
}  // namespace dca

#endif  // TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_SINGLE_SECTOR_ACCUMULATION_TEST_HPP
