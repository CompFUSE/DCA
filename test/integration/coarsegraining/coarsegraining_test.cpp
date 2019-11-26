// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Run the coarsegraining on a two-band model with n-sites, and compare the results with the
// precomputed result from a singleband model with 2*n sites representing the same physical system.
// The Green's function for the two systems should be equal in real space, and related by the
// transformation G_singleband(k) = G_aa(k) + G_ab(k) and  G_singleband(k+\pi) = G_aa(k) - G_ab(k)
// in momentum space, where G_aa is the propagator between the same orbitals, and G_ab between the
// two different orbitals.
//
// Note: as the code does not support reinitialization of the cluster with a different model or
//       size, the singleband results are stored as a baseline in an hdf5 file, while the
//       two-band model is run during this test. If the baseline needs to be update, please define
//       the flag UPDATE_BASELINE and run the resulting executable.
#undef UPDATE_BASELINE

#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_sp.hpp"

#include <iostream>
#include <string>

#include "dca/config/cmake_options.hpp"
#include "dca/config/threading.hpp"
//
#include "gtest/gtest.h"
//
#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/singleband_chain.hpp"
#include "dca/phys/models/analytic_hamiltonians/twoband_chain.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"

// Set to true to dump the result in an hdf5 file.
constexpr bool write_G_r_w = false;

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/coarsegraining/";

using Concurrency = dca::parallel::MPIConcurrency;
using Model1 = dca::phys::models::TightBindingModel<
    dca::phys::models::singleband_chain<dca::phys::domains::no_symmetry<2>>>;
using Model2 = dca::phys::models::TightBindingModel<
    dca::phys::models::twoband_chain<dca::phys::domains::no_symmetry<2>>>;

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
using Coarsegraining = dca::phys::clustermapping::CoarsegrainingSp<Parameters>;

using RDmn = Data::RClusterDmn;
using KDmn = Data::KClusterDmn;

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
    cluster_mapping_obj.compute_G_K_w(data.Sigma_lattice, data.G_k_w);
  }
  else {
    computeMockSigma(data.Sigma);
    cluster_mapping_obj.compute_G_K_w(data.Sigma, data.G_k_w);
  }

  dca::math::transform::FunctionTransform<KDmn, RDmn>::execute(data.G_k_w, data.G_r_w);

  if (concurrency.id() == 0) {
    const std::string baseline_name = test_dca_plus ? "coarsegraining_dca_plus_baseline.hdf5"
                                                    : "coarsegraining_dca_baseline.hdf5";
    std::vector<std::complex<double>> raw_data;

#ifdef UPDATE_BASELINE
    raw_data.resize(data.G_r_w.size());
    std::copy_n(data.G_r_w.values(), data.G_r_w.size(), raw_data.data());

    dca::io::HDF5Writer writer;
    writer.open_file(input_dir + baseline_name);
    writer.execute(data.G_r_w.get_name(), raw_data);
    writer.close_file();

    if (write_G_r_w && test_dca_plus == false) {
      writer.open_file("coarse_singleband.hdf5");
      data.write(writer);
      parameters.write(writer);
    }
#else
    using SingleBandBDmn = dca::func::dmn_0<dca::func::dmn<1>>;
    using SDmn = dca::func::dmn_0<dca::func::dmn<2>>;
    using WDmn = dca::func::dmn_0<dca::phys::domains::frequency_domain>;
    using SingleBandRDmn = dca::func::dmn_0<dca::func::dmn<6>>;
    using SinglebandNuDmn = dca::func::dmn_variadic<SingleBandBDmn, SDmn>;
    dca::func::function<std::complex<double>,
                        dca::func::dmn_variadic<SinglebandNuDmn, SinglebandNuDmn, SingleBandRDmn, WDmn>>
        G_check(data.G_r_w.get_name());

    dca::io::HDF5Reader reader;
    reader.open_file(input_dir + baseline_name);

    reader.execute(G_check.get_name(), raw_data);
    reader.close_file();
    ASSERT_EQ(raw_data.size(), G_check.size());
    std::copy_n(raw_data.data(), raw_data.size(), G_check.values());

    if (write_G_r_w && test_dca_plus == false) {
      dca::io::HDF5Writer writer;
      writer.open_file("coarse_bilayer.hdf5");
      data.write(writer);
      parameters.write(writer);
    }

    // The Greens function in real space should be the same up to a relabelling.
    for (int w = 0; w < WDmn::dmn_size(); ++w)
      for (int r = 0; r < RDmn::dmn_size(); ++r)
        for (int s = 0; s < SDmn::dmn_size(); ++s) {
          EXPECT_LE(std::abs(data.G_r_w(r % 2, s, 0, s, r / 2, w) - G_check(0, s, 0, s, r, w)), 1e-7);
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
