// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Compare the coarsegraining of a bileyer model with one site, with a singleband model with two
// sites.

#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_sp.hpp"

#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/singleband_chain.hpp"
#include "dca/phys/models/analytic_hamiltonians/twoband_chain.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"

// If this flag is defined, prepare the baseline using the singleband model, otherwise run the test
// comparing the bilayer model with the baseline.
#undef UPDATE_BASELINE

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/coarsegraining/";

using Concurrency = dca::parallel::MPIConcurrency;
using Model1 = dca::phys::models::TightBindingModel<
    dca::phys::models::singleband_chain<dca::phys::domains::no_symmetry<2>>>;
using Model2 = dca::phys::models::TightBindingModel<
    dca::phys::models::twoband_chain<dca::phys::domains::no_symmetry<2>>>;
using Threading = dca::parallel::stdthread;

#ifdef UPDATE_BASELINE
using Parameters = dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler,
                                                 Model1, void, dca::phys::solver::CT_AUX>;
const std::string input = input_dir + "input_singleband.json";
#else
using Parameters = dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler,
                                                 Model2, void, dca::phys::solver::CT_AUX>;
const std::string input = input_dir + "input_bilayer.json";
#endif  // UPDATE_BASELINE

using Data = dca::phys::DcaData<Parameters>;
using Coarsegraining =
    dca::phys::clustermapping::CoarsegrainingSp<Parameters>;

template <class SigmaType>
void computeMockSigma(SigmaType& Sigma);

void performTest(const bool test_dca_plus) {
  static Concurrency concurrency(0, nullptr);

  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input);

  static bool model_initialized = false;
  if (!model_initialized) {
    parameters.update_model();
    parameters.update_domains();
    model_initialized = true;
  }

  Data data(parameters);
  data.initialize();

  Coarsegraining cluster_mapping_obj(parameters);

  if (test_dca_plus) {
    computeMockSigma(data.Sigma_lattice);
    cluster_mapping_obj.compute_G_K_w(data.H_HOST, data.Sigma_lattice, data.G_k_w);
  }
  else {
    computeMockSigma(data.Sigma);
    cluster_mapping_obj.compute_G_K_w(data.H_HOST, data.Sigma, data.G_k_w);
  }

  if (concurrency.id() == 0) {
    const std::string baseline_name = test_dca_plus ? "coarsegraining_dca_plus_baseline.hdf5"
                                                    : "coarsegraining_dca_baseline.hdf5";
    std::vector<std::complex<double>> raw_data;

#ifdef UPDATE_BASELINE
    raw_data.resize(data.G_k_w.size());
    std::copy_n(data.G_k_w.values(), data.G_k_w.size(), raw_data.data());

    dca::io::HDF5Writer writer;
    writer.open_file(input_dir + baseline_name);
    writer.execute(data.G_k_w.get_name(), raw_data);
    writer.close_file();
#else
    using SingleBandBDmn = dca::func::dmn_0<dca::func::dmn<1>>;
    using SDmn = dca::func::dmn_0<dca::func::dmn<2>>;
    using WDmn = dca::func::dmn_0<dca::phys::domains::frequency_domain>;
    using SingleBandKDmn = dca::func::dmn_0<dca::func::dmn<2>>;
    using SinglebandNuDmn = dca::func::dmn_variadic<SingleBandBDmn, SDmn>;
    dca::func::function<std::complex<double>,
                        dca::func::dmn_variadic<SinglebandNuDmn, SinglebandNuDmn, SingleBandKDmn, WDmn>>
        G_check(data.G_k_w.get_name());

    dca::io::HDF5Reader reader;
    reader.open_file(input_dir + baseline_name);

    reader.execute(G_check.get_name(), raw_data);
    reader.close_file();
    ASSERT_EQ(raw_data.size(), G_check.size());
    std::copy_n(raw_data.data(), raw_data.size(), G_check.values());

    // The desired result is G(k=0) = G_aa(k=0) + G_ab(k=0) and G(k=\pi) = G_aa(k=0) - G_ab(k=0).
    // This follows from c(0) = (c_aa(0) + c_ab(0)) / sqrt(2)
    // and c(\pi) = (c_aa(0) - c_ab(0)) / sqrt(2).
    for (int w = 0; w < WDmn::dmn_size(); ++w)
      for (int s = 0; s < SDmn::dmn_size(); ++s) {
        const std::complex<double> G_k_0 = G_check(0, s, 0, s, 0, w);
        const std::complex<double> G_k_pi = G_check(0, s, 0, s, 1, w);
        const std::complex<double> G_aa = (G_k_0 + G_k_pi) / 2.;
        const std::complex<double> G_ab = (G_k_0 - G_k_pi) / 2.;
        EXPECT_LE(std::abs(G_aa - data.G_k_w(0, s, 0, s, 0, w)), 1e-7);
        EXPECT_LE(std::abs(G_aa - data.G_k_w(1, s, 1, s, 0, w)), 1e-7);
        EXPECT_LE(std::abs(G_ab - data.G_k_w(0, s, 1, s, 0, w)), 1e-7);
        EXPECT_LE(std::abs(G_ab - data.G_k_w(1, s, 0, s, 0, w)), 1e-7);
      }
#endif  // UPDATE_BASELINE
  }
}

TEST(CoarsegrainingTest, DCABilayerVsSingleband) {
  performTest(false);
}

TEST(CoarsegrainingTest, DCAPlusBilayerVsSingleband) {
  performTest(true);
}

template <class SigmaType>
void computeMockSigma(SigmaType& Sigma) {
  using BDmn = dca::phys::domains::electron_band_domain;
  using WDmn = dca::phys::domains::frequency_domain;
  const int n_k = Sigma[4];

  const double U = 4.;
  const std::complex<double> imag(0, 1.);

  for (int w = 0; w < WDmn::get_size(); ++w) {
    const double w_val = WDmn::get_elements()[w];
    const std::complex<double> sigma_val = U * U / (4. * imag * w_val);
    for (int k = 0; k < n_k; ++k)
      for (int s = 0; s < 2; ++s)
        for (int b = 0; b < BDmn::get_size(); ++b)
          Sigma(b, s, b, s, k, w) = sigma_val;
  }
}
