// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class contains all functions needed for the MOMS DCA calculation.
//
// TODO: Interpolate Sigma if \beta from file != \beta ?

#ifndef DCA_PHYS_DCA_DATA_DCA_DATA_HPP
#define DCA_PHYS_DCA_DATA_DCA_DATA_HPP

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstring>
#include <iostream>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>
#include <variant>

#include "dca/distribution/dist_types.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/util/real_complex_conversion.hpp"
#include "dca/io/reader.hpp"
#include "dca/io/writer.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/dca_algorithms/compute_band_structure.hpp"
#include "dca/phys/dca_algorithms/compute_free_greens_function.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/domains/cluster/interpolation/hspline_interpolation/hspline_interpolation.hpp"
#include "dca/phys/domains/cluster/momentum_exchange_domain.hpp"
#include "dca/phys/domains/quantum/brillouin_zone_cut_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_exchange_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/error_computation_type.hpp"
#include "dca/phys/four_point_type.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/models/traits.hpp"
#include "dca/util/timer.hpp"
#include "dca/util/to_string.hpp"
#include "dca/distribution/dist_types.hpp"
#ifdef DCA_WITH_ADIOS2
#include "dca/io/adios2/adios2_writer.hpp"
#endif
namespace dca {
namespace phys {
// dca::phys::

template <class Parameters, DistType DIST = DistType::NONE>
class DcaData {
public:
  static constexpr DistType DT = DIST;
  using profiler_type = typename Parameters::profiler_type;

  using Concurrency = typename Parameters::concurrency_type;
  using Lattice = typename Parameters::lattice_type;
  constexpr static int DIMENSION = Lattice::DIMENSION;

  using Real = typename Parameters::Real;
  using Scalar = typename Parameters::Scalar;
  using TpAccumulatorPrec = typename Parameters::TPAccumPrec;
  using TpComplex = std::complex<TpAccumulatorPrec>;
  using TDmn = func::dmn_0<domains::time_domain>;
  using WDmn = func::dmn_0<domains::frequency_domain>;
  using WVertexDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using WExchangeDmn = func::dmn_0<domains::FrequencyExchangeDomain>;

  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;  // orbital-spin index
  using NuNuDmn = func::dmn_variadic<NuDmn, NuDmn>;

  using CDA = ClusterDomainAliases<DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterType = typename CDA::KClusterType;
  using KClusterDmn = typename CDA::KClusterDmn;
  using RHostDmn = typename CDA::RSpHostDmn;
  using KHostDmn = typename CDA::KSpHostDmn;
  using KExchangeDmn = func::dmn_0<domains::MomentumExchangeDomain>;

  using KCutDmn = func::dmn_0<domains::brillouin_zone_cut_domain<101>>;

  using NuKCutDmn = func::dmn_variadic<NuDmn, KCutDmn>;
  using NuNuKWDmn = func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>;

  using SpGreensFunction =
      func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>>;
  using SpRGreensFunction =
      func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>;

  using TpGreensFunction =
      func::function<TpComplex,
                     func::dmn_variadic<BDmn, BDmn, BDmn, BDmn, KClusterDmn, WVertexDmn,
                                        KClusterDmn, WVertexDmn, KExchangeDmn, WExchangeDmn>,
                     DT>;

  DcaData(Parameters& parameters_ref);

  /** These reads are used by analysis programs only for now.
   */
  void read(const std::string& filename);

  /** prefer this as it allows for more sensible handling of appendable files like bp4.
   *  currently this is used by chi_q_omega only.
   */
  void read(dca::io::Reader<Concurrency>& reader);

  template <typename Writer>
  void write(Writer& writer);

#ifdef DCA_WITH_ADIOS2
  void writeDistributedG4Adios(io::ADIOS2Writer<Concurrency>& writer);
#endif

  void initialize();
  void initializeH0_and_H_i();
  void initialize_G0();
  void initializeSigma(const std::string& filename);
  void readSigmaFile(io::Reader<Concurrency>& reader);

  void compute_single_particle_properties();
  void compute_Sigma_bands();

  void print_Sigma_QMC_versus_Sigma_cg();

private:
  Parameters& parameters_;
  Concurrency& concurrency_;

public:
  func::function<int, NuNuDmn> H_symmetry;

  // Interaction Hamiltonian. Each entry H_interactions(nu1, nu2, delta_r) represents the
  // correlation strength between the two orbitals nu1 nu2 at distance delta_r. This correlation
  // must be symmetric, or double counted, i.e.
  // H_interactions(nu1, nu2, delta_r) == H_interactions(nu2, nu1, -delta_r). Each pair of terms
  // represents a single addendum in the physical hamiltonian proportional to n_{nu1} * n_{nu2}, or
  // H = \sum_{nu1, nu2, r1, r2} H_interactions(nu1, nu2, r1 - r2) n_{nu1} n_{nu2} / 2.
  func::function<Real, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn>> H_interactions;

  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn>> H_DCA;
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn>> H_HOST;

  func::function<Real, NuKCutDmn> band_structure;

  func::function<std::complex<Real>, NuKCutDmn> Sigma_band_structure;

  func::function<std::complex<Real>, NuKCutDmn> Sigma_cluster_band_structure;
  func::function<std::complex<Real>, NuKCutDmn> Sigma_lattice_band_structure;

  func::function<std::complex<Real>, NuKCutDmn> Sigma_band_structure_interpolated;
  func::function<std::complex<Real>, NuKCutDmn> Sigma_band_structure_coarsegrained;

  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn>> G_k;  //("Greens-k-lattice");
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn>> S_k;  //("Sigma-k-lattice");
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, RHostDmn>> S_r;  //("Sigma-r-lattice");

  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>> Sigma;

  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>> Sigma_cluster;
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn, WDmn>> Sigma_lattice;
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn, WDmn>>
      Sigma_lattice_interpolated;
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn, WDmn>>
      Sigma_lattice_coarsegrained;

  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>> G_k_w;
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, TDmn>> G_k_t;
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>> G_r_w;
  func::function<Scalar, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>> G_r_t;

  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>> G0_k_w;
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, TDmn>> G0_k_t;
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>> G0_r_w;
  func::function<Scalar, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>> G0_r_t;

  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>>
      G0_k_w_cluster_excluded;
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, TDmn>>
      G0_k_t_cluster_excluded;
  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>
      G0_r_w_cluster_excluded;
  func::function<Scalar, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>> G0_r_t_cluster_excluded;

  func::function<Real, NuDmn> orbital_occupancy;

public:  // Optional members getters.
  auto& get_G_k_w_error() {
    if (not G_k_w_err_)
      G_k_w_err_.reset(new SpGreensFunction("G_k_w-error"));
    return *G_k_w_err_;
  }
  auto& get_G_r_w_error() {
    if (not G_r_w_err_)
      G_r_w_err_ = std::make_unique<SpRGreensFunction>("G_r_w-error");
    return *G_r_w_err_;
  }
  auto& get_G_k_w_stdv() {
    if (not G_k_w_err_)
      G_k_w_err_.reset(new SpGreensFunction("cluster_greens_function_G_k_w-stddev"));
    return *G_k_w_err_;
  }
  auto& get_Sigma_stdv() {
    if (not Sigma_err_)
      Sigma_err_.reset(new SpGreensFunction("Self_Energy-stddev"));
    return *Sigma_err_;
  }
  auto& get_Sigma_error() {
    if (not Sigma_err_)
      Sigma_err_.reset(new SpGreensFunction("Self_Energy-error"));
    return *Sigma_err_;
  }
  auto& get_G4() {
    assert(!G4_.empty());
    return G4_;
  }
  auto& get_G4_error() {
    assert(!G4_err_.empty());
    return G4_err_;
  }
  auto& get_G4_stdv() {
    assert(!G4_err_.empty());
    return G4_err_;
  }

  // The non density-density Hamiltonian is given by:
  // H = \sum(nu1, nu2, nu3, nu4, r1, r2) c^+(nu1, r1) c(nu2, r1) c^+(nu3, r2) c(nu4, r2) *
  //     non_density_interactions_(nu1, nu2, nu3, nu4, r1 - r2)
  // Note: this contribution to the Hamiltonian is not Real counted.
  auto& get_non_density_interactions() {
    if (not non_density_interactions_) {
      non_density_interactions_ = std::make_unique<
          func::function<Real, func::dmn_variadic<NuDmn, NuDmn, NuDmn, NuDmn, RClusterDmn>>>(
          "non_density_interaction");
    }
    return *non_density_interactions_;
  }
  const auto& get_non_density_interactions() const {
    assert(non_density_interactions_);
    return *non_density_interactions_;
  }

  bool has_non_density_interactions() const {
    return (bool)non_density_interactions_;
  }

private:  // Optional members.
  std::unique_ptr<SpGreensFunction> G_k_w_err_;
  std::unique_ptr<SpRGreensFunction> G_r_w_err_;
  std::unique_ptr<SpGreensFunction> Sigma_err_;
  std::vector<TpGreensFunction> G4_;
  std::vector<TpGreensFunction> G4_err_;
  std::unique_ptr<func::function<Real, func::dmn_variadic<NuDmn, NuDmn, NuDmn, NuDmn, RClusterDmn>>>
      non_density_interactions_;
};

template <class Parameters, DistType DT>
DcaData<Parameters, DT>::DcaData(/*const*/ Parameters& parameters_ref)
    : parameters_(parameters_ref),
      concurrency_(parameters_.get_concurrency()),

      H_symmetry("H_symmetry"),
      H_interactions("interaction-matrix"),

      H_DCA("H_DCA"),
      H_HOST("H_HOST"),

      band_structure("band-structure"),

      Sigma_band_structure("Sigma-band-structure"),

      Sigma_cluster_band_structure("Sigma-cluster-band-structure"),
      Sigma_lattice_band_structure("Sigma-lattice-band-structure"),

      Sigma_band_structure_interpolated("Sigma-band-structure-interpolated"),
      Sigma_band_structure_coarsegrained("Sigma-band-structure-coarsegrained"),

      G_k("Greens-k-lattice"),
      S_k("Sigma-k-lattice"),
      S_r("Sigma-r-lattice"),

      Sigma("Self_Energy"),

      Sigma_cluster("Self-Energy-cluster"),
      Sigma_lattice("Self-energy-lattice"),

      Sigma_lattice_interpolated("Sigma_lattice_interpolated"),
      Sigma_lattice_coarsegrained("Sigma_lattice_coarsegrained"),

      G_k_w("cluster_greens_function_G_k_w"),
      G_k_t("cluster_greens_function_G_k_t"),
      G_r_w("cluster_greens_function_G_r_w"),
      G_r_t("cluster_greens_function_G_r_t"),

      G0_k_w("free_cluster_greens_function_G0_k_w"),
      G0_k_t("free_cluster_greens_function_G0_k_t"),
      G0_r_w("free_cluster_greens_function_G0_r_w"),
      G0_r_t("free_cluster_greens_function_G0_r_t"),

      G0_k_w_cluster_excluded("cluster_excluded_greens_function_G0_k_w"),
      G0_k_t_cluster_excluded("cluster_excluded_greens_function_G0_k_t"),
      G0_r_w_cluster_excluded("cluster_excluded_greens_function_G0_r_w"),
      G0_r_t_cluster_excluded("cluster_excluded_greens_function_G0_r_t"),

      orbital_occupancy("orbital_occupancy") {
  H_symmetry = -1;

  // Reserve storage in advance such that we don't have to copy elements when we fill the vector.
  // We want to avoid copies because function's copy ctor does not copy the name (and because copies
  // are expensive).
  for (auto channel : parameters_.get_four_point_channels()) {
    // Allocate memory for G4, eventually distributed among all processes.
    if (parameters_.get_g4_distribution() == DistType::BLOCKED) {
      G4_.emplace_back("G4_" + toString(channel), concurrency_);
      G4_err_.emplace_back("G4_" + toString(channel) + "_err", concurrency_);
    }
    else {
      G4_.emplace_back("G4_" + toString(channel));
      G4_err_.emplace_back("G4_" + toString(channel) + "_err");
    }
  }
}

template <class Parameters, DistType DT>
void DcaData<Parameters, DT>::read(const std::string& filename) {
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "\n\n\t starts reading \n\n";
    dca::io::Reader<typename Parameters::concurrency_type> reader(concurrency_,
                                                                  parameters_.get_output_format());
    reader.open_file(filename);
    read(reader);
    reader.close_file();
  }

  concurrency_.broadcast(parameters_.get_chemical_potential());
#ifndef NDEBUG
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "broadcasted chemical potential: " << parameters_.get_chemical_potential() << '\n';
  }
#endif
  concurrency_.broadcast_object(Sigma);
#ifndef NDEBUG
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "broadcasted Sigma \n";
  }
#endif

  if (parameters_.isAccumulatingG4()) {
    concurrency_.broadcast_object(G_k_w);
#ifndef NDEBUG
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "broadcasted G_k_w \n";
  }
#endif

  for (auto& G4_channel : G4_) {
      concurrency_.broadcast_object(G4_channel);
#ifndef NDEBUG
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "broadcasted G4_channel \n";
  }
#endif
  }
  }
}

template <class Parameters, DistType DT>
void DcaData<Parameters, DT>::read(dca::io::Reader<typename Parameters::concurrency_type>& reader) {
  std::size_t step_count = reader.getStepCount();
  for (std::size_t i = 0; i < step_count; ++i) {
    reader.begin_step();
  }

  reader.open_group("parameters");

  reader.open_group("physics");
  reader.execute("chemical-potential", parameters_.get_chemical_potential());
  reader.close_group();

  reader.close_group();

  reader.open_group("functions");

  reader.execute(Sigma);

  if (parameters_.isAccumulatingG4()) {
    std::cout << "Trying to read Gkw since we are accumulating G4\n";
    reader.execute(G_k_w);

    // Try to read G4 with a legacy name.
    if (parameters_.get_four_point_channels().size() == 1) {
      reader.execute("G4", G4_[0]);
    }

    for (auto& G4_channel : G4_)
      reader.execute(G4_channel);
  }

  reader.close_group();
}

#ifdef DCA_WITH_ADIOS2
template <class Parameters, DistType DIST>
void DcaData<Parameters, DIST>::writeDistributedG4Adios(io::ADIOS2Writer<Concurrency>& writer) {
  if constexpr (DIST == DistType::BLOCKED || DIST == DistType::LINEAR) {
    if (parameters_.isAccumulatingG4() && parameters_.get_g4_output_format() == "ADIOS2" &&
        parameters_.get_g4_distribution() != DistType::NONE) {
      std::cerr << "trying to write G4 to adios on rank: " << concurrency_.id() << '\n';

      for (const auto& G4_channel : G4_) {
#ifndef NDEBUG
        std::cerr << "Writing G4_channel:" << G4_channel.get_name()
                  << "on rank: " << concurrency_.id() << '\n';
        std::cerr << "start: " << G4_channel.get_start() << "   end: " << G4_channel.get_end()
                  << '\n';
        auto str_sub_ind_start = vectorToString(G4_channel.get_start_subindex());
        auto str_sub_ind_end = vectorToString(G4_channel.get_end_subindex());
        std::cerr << "start subind: " << str_sub_ind_start << "   end: " << str_sub_ind_end << '\n';
#endif
        writer.execute(G4_channel);
      }
    }
  }
}
#endif

template <class Parameters, DistType DT>
template <typename Writer>
void DcaData<Parameters, DT>::write(Writer& writer) {
  writer.open_group("functions");

  writer.execute(band_structure);

  if (parameters_.do_dca_plus()) {
    writer.execute(Sigma_band_structure);
    writer.execute(Sigma_cluster_band_structure);
    writer.execute(Sigma_lattice_band_structure);
    writer.execute(Sigma_band_structure_interpolated);
    writer.execute(Sigma_band_structure_coarsegrained);

    writer.execute(S_k);
    writer.execute(S_r);

    writer.execute(G_k);
  }

  else {
    // Compute Sigma-r-DCA for the lowest frequency via Fourier transformation of DCA cluster Sigma.
    func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn>> S_r_DCA(
        "Sigma-r-DCA");

    func::function<std::complex<Real>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn>> S_k_DCA(
        "Sigma-k-DCA");
    std::memcpy(
        &S_k_DCA(0), &Sigma(0, 0, 0, WDmn::dmn_size() / 2),
        sizeof(std::complex<Real>) * std::pow(2 * BDmn::dmn_size(), 2.) * KClusterDmn::dmn_size());
    math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(S_k_DCA, S_r_DCA);

    writer.execute(S_r_DCA);
  }

  writer.execute(Sigma);
  writer.execute(Sigma_err_);
  writer.execute(Sigma_cluster);

  if (parameters_.dump_lattice_self_energy()) {
    if (parameters_.do_dca_plus())
      writer.execute(Sigma_lattice);
    else if (parameters_.doPostInterpolation())
      writer.execute(Sigma_lattice_interpolated);
  }

  if (parameters_.dump_cluster_Greens_functions()) {
    writer.execute(G_k_w);
    writer.execute(G_k_w_err_);
    writer.execute(G_r_w);
    writer.execute(G_r_w_err_);
    writer.execute(G_k_t);
    writer.execute(G_r_t);

    writer.execute(G0_k_w);
    writer.execute(G0_r_w);
    writer.execute(G0_k_t);
    writer.execute(G0_r_t);

    writer.execute(G0_k_w_cluster_excluded);
    writer.execute(G0_r_w_cluster_excluded);
    writer.execute(G0_k_t_cluster_excluded);
    writer.execute(G0_r_t_cluster_excluded);
  }
  else {
    writer.execute(G_k_w);
    writer.execute(G_k_w_err_);
  }

  // When distributed_g4_enabled, one should assume G4 size is fairly large and then should not
  // accumulate G4 into one node and thus cannot write it out
  if (parameters_.isAccumulatingG4()) {
    if constexpr (DT != DistType::BLOCKED) {
      for (const auto& G4_channel : G4_)
        writer.execute(G4_channel);

      if (parameters_.get_error_computation_type() != ErrorComputationType::NONE) {
        for (const auto& G4_channel_err : G4_err_)
          writer.execute(G4_channel_err);
      }
    }
#ifdef DCA_WITH_ADIOS2
    else {
      // special adios writer only for block or linear distributed G4
      if (writer.isADIOS2())
        writeDistributedG4Adios(std::get<io::ADIOS2Writer<Concurrency>>(writer.getUnderlying()));
    }
#endif
  }
  writer.close_group();
}

template <class Parameters, DistType DT>
void DcaData<Parameters, DT>::initialize() {
  initializeH0_and_H_i();
  initialize_G0();
}

template <class Parameters, DistType DT>
void DcaData<Parameters, DT>::initializeH0_and_H_i() {
  util::Timer("H_0 and H_int initialization", concurrency_.id() == concurrency_.first());

  Parameters::model_type::initializeH0(parameters_, H_DCA);
  Parameters::model_type::initializeH0(parameters_, H_HOST);

  Parameters::model_type::initializeHInteraction(H_interactions, parameters_);

  // Check symmetry of H_interactions.
  const int r0 = RClusterDmn::parameter_type::origin_index();
  for (int r = 0; r < RClusterDmn::dmn_size(); ++r) {
    const int minus_r = RClusterDmn::parameter_type::subtract(r, r0);
    for (int nu2 = 0; nu2 < NuDmn::dmn_size(); ++nu2)
      for (int nu1 = 0; nu1 < NuDmn::dmn_size(); ++nu1) {
        if (std::abs(H_interactions(nu1, nu2, r) - H_interactions(nu2, nu1, minus_r)) > 1e-8) {
          std::cout << r << " , " << minus_r << " , " << H_interactions(nu1, nu2, r) << " , " << H_interactions(nu2, nu1, minus_r) << "\n";
          throw(std::logic_error("Double counting is not consistent."));
        }
      }
  }

  if constexpr (models::HasInitializeNonDensityInteractionMethod<Parameters>::value) {
    models::initializeNonDensityInteraction<Parameters>(get_non_density_interactions(), parameters_);
  }

  Parameters::model_type::initialize_H_symmetries(H_symmetry);

  compute_band_structure<Parameters>::execute(parameters_, band_structure);
}

template <class Parameters, DistType DT>
void DcaData<Parameters, DT>::initialize_G0() {
  try {
    profiler_type prof(__FUNCTION__, "DcaData", __LINE__);

    util::Timer("G_0 initialization", concurrency_.id() == concurrency_.first());

    // Compute G0_k_w.
    compute_G0_k_w(H_DCA, parameters_.get_chemical_potential(),
                   parameters_.get_coarsegraining_threads(), G0_k_w);
    Symmetrize<Parameters>::execute(G0_k_w, H_symmetry, true);

    // Compute G0_k_t.
    compute_G0_k_t(H_DCA, parameters_.get_chemical_potential(), parameters_.get_beta(), G0_k_t);
    Symmetrize<Parameters>::execute(G0_k_t, H_symmetry, true);

    // Compute G0_r_w.
    math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(G0_k_w, G0_r_w);
    Symmetrize<Parameters>::execute(G0_r_w, H_symmetry, true);

    // Compute G0_r_t.
    math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(G0_k_t, G0_r_t);
    Symmetrize<Parameters>::execute(G0_r_t, H_symmetry, true);

    // Initialize the cluster excluded Green's functions with the corresponding free Green's
    // functions.
    G0_k_w_cluster_excluded = G0_k_w;
    G0_k_t_cluster_excluded = G0_k_t;
    G0_r_w_cluster_excluded = G0_r_w;
    G0_r_t_cluster_excluded = G0_r_t;
  }
  catch (const std::exception& exc) {
    std::throw_with_nested(std::runtime_error("Failure in initialization of G0!"));
  }
}

template <class Parameters, DistType DT>
void DcaData<Parameters, DT>::initializeSigma(const std::string& filename) {
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "reading Sigma File\n";
    io::IOType sigma_file_io = io::extensionToIOType(filename);
    io::Reader reader(concurrency_, sigma_file_io);
    if (sigma_file_io == io::IOType::HDF5 || sigma_file_io == io::IOType::ADIOS2) {
      io::Reader reader(concurrency_, sigma_file_io);
      int hdf5_last_iteration = -1;
      reader.open_file(filename);
      long step_count = reader.getStepCount();
      // Work around odd way hdf5 steps get written
      int completed_iteration = 0;
      std::cout << step_count << " steps found.\n";
      if (step_count >= 0) {
      find_step:
        for (std::size_t i = 0; i < step_count; ++i) {
          reader.begin_step();
          std::cerr << "current step " << i << '\n';
          bool has_iteration =
              reader.execute("DCA-loop-functions/completed-iteration", completed_iteration);
          std::cerr << "completed_iteration " << completed_iteration << '\n';
          if (has_iteration && (i > completed_iteration)) {
            std::cerr << "past complete iterations " << completed_iteration << "at step " << i
                      << '\n';
            hdf5_last_iteration = i;
            reader.close_file();
            reader.open_file(filename);
            step_count = hdf5_last_iteration;
            goto find_step;
          }
          if (i < step_count - 1)
            reader.end_step();
        }
        readSigmaFile(reader);
        reader.end_step();
      }
      else {
        readSigmaFile(reader);
      }
    }
    else {
      reader.open_file(filename);
      std::size_t step_count = reader.getStepCount();
      for (std::size_t i = 0; i < step_count; ++i) {
        reader.begin_step();
        reader.end_step();
      }
      readSigmaFile(reader);
      reader.close_file();
    }
  }
  concurrency_.broadcast(parameters_.get_chemical_potential());
  #ifndef NDEBUG
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "broadcasted chemical potential: " << parameters_.get_chemical_potential();
  }
#endif
  concurrency_.broadcast(Sigma);
  #ifndef NDEBUG
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "broadcasted Sigma \n";
  }
#endif
}

template <class Parameters, DistType DT>
void DcaData<Parameters, DT>::readSigmaFile(io::Reader<Concurrency>& reader) {
  reader.open_group("DCA-loop-functions");
  std::vector<Real> chemical_potentials;
  bool chemical_potential_present = reader.execute("chemical-potential", chemical_potentials);
  std::vector<int> completed_iterations;
  bool has_iteration = reader.execute("completed-iteration", completed_iterations);
  if (chemical_potential_present && has_iteration) {
    int completed_iteration =
        *std::max_element(completed_iterations.begin(), completed_iterations.end());
    std::cout << "chemical-potential from Sigma file: " << chemical_potentials[completed_iteration]
              << '\n';
    parameters_.get_chemical_potential() = chemical_potentials[completed_iteration];
  }
  else {
    throw std::runtime_error(
        "readSigmaFile failed, initial-self-energy file is missing chemical potential data!");
  }
  reader.close_group();

  reader.open_group("functions");
  reader.execute(Sigma);
  reader.close_group();
}

template <class Parameters, DistType DT>
void DcaData<Parameters, DT>::compute_single_particle_properties() {
  {
    std::memcpy(
        &S_k(0), &Sigma_lattice(0, 0, 0, WDmn::dmn_size() / 2),
        sizeof(std::complex<Real>) * std::pow(2 * BDmn::dmn_size(), 2.) * KHostDmn::dmn_size());

    math::transform::FunctionTransform<KHostDmn, RHostDmn>::execute(S_k, S_r);
  }

  {
    int w_ind = WDmn::dmn_size() / 2;

    dca::linalg::Matrix<std::complex<Real>, dca::linalg::CPU> I_k("I_matrix", NuDmn::dmn_size());
    dca::linalg::Matrix<std::complex<Real>, dca::linalg::CPU> G_inv("G_inv", NuDmn::dmn_size());

    // Allocate the work space for inverse only once.
    dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
    dca::linalg::Vector<std::complex<Real>, dca::linalg::CPU> work;

    std::complex<Real> i_wm_plus_mu;

    i_wm_plus_mu.real(parameters_.get_chemical_potential());
    i_wm_plus_mu.imag(WDmn::get_elements()[w_ind]);

    for (int i = 0; i < NuDmn::dmn_size(); i++)
      I_k(i, i) = i_wm_plus_mu;

    for (int k_ind = 0; k_ind < KHostDmn::dmn_size(); k_ind++) {
      for (int j = 0; j < NuDmn::dmn_size(); j++)
        for (int i = 0; i < NuDmn::dmn_size(); i++)
          G_inv(i, j) = I_k(i, j) - H_HOST(i, j, k_ind) - Sigma_lattice(i, j, k_ind, w_ind);

      dca::linalg::matrixop::inverse(G_inv, ipiv, work);

      for (int j = 0; j < NuDmn::dmn_size(); j++)
        for (int i = 0; i < NuDmn::dmn_size(); i++)
          G_k(i, j, k_ind) = G_inv(i, j);
    }
  }
}

template <class Parameters, DistType DT>
void DcaData<Parameters, DT>::compute_Sigma_bands() {
  {
    // We now always sum over both spin sectors so this is appropriate for both real and complex g0
    Sigma_band_structure.reset();
    Sigma_cluster_band_structure.reset();

    std::vector<std::pair<Real, int>> length_and_distance(KClusterDmn::dmn_size(),
                                                          std::pair<Real, int>(0, -1));

    auto kConvert = [](auto& kvec) -> std::vector<Real> {
      std::vector<Real> k_converted(kvec.size());
      std::transform(kvec.begin(), kvec.end(), k_converted.begin(),
                     [](auto& val) -> typename decltype(k_converted)::value_type {
                       return static_cast<typename decltype(k_converted)::value_type>(val);
                     });
      return k_converted;
    };

    if constexpr (!std::is_same_v<typename KClusterType::Scalar, Real>) {
      auto full_prec_super_basis = KClusterType::get_super_basis_vectors();
      std::vector<std::vector<Real>> ksuper_basis(full_prec_super_basis.size());
      for (int i = 0; i < full_prec_super_basis.size(); ++i) {
        ksuper_basis[i] = kConvert(full_prec_super_basis[i]);
      }
      for (int k_ind = 0; k_ind < KCutDmn::dmn_size(); ++k_ind) {
        auto kcut_slice = kConvert(KCutDmn::get_elements()[k_ind]);
        auto k_vec = domains::cluster_operations::translate_inside_cluster(kcut_slice, ksuper_basis);

        for (int K_ind = 0; K_ind < KClusterDmn::dmn_size(); ++K_ind) {
          auto kdmn_slice = kConvert(KClusterDmn::get_elements()[K_ind]);
          length_and_distance[K_ind].second = K_ind;
          length_and_distance[K_ind].first =
              domains::cluster_operations::minimal_distance(k_vec, kdmn_slice, ksuper_basis);
        }
        std::sort(length_and_distance.begin(), length_and_distance.end());

        int result_ind = length_and_distance[0].second;

        for (int nu_ind = 0; nu_ind < 2 * BDmn::dmn_size(); ++nu_ind) {
          Sigma_band_structure(nu_ind, k_ind) =
              Sigma(nu_ind, nu_ind, result_ind, WDmn::dmn_size() / 2);
          Sigma_cluster_band_structure(nu_ind, k_ind) =
              Sigma_cluster(nu_ind, nu_ind, result_ind, WDmn::dmn_size() / 2);
        }
      }
    }
    else {
      for (int k_ind = 0; k_ind < KCutDmn::dmn_size(); ++k_ind) {
        std::vector<double> k_vec = domains::cluster_operations::translate_inside_cluster(
            KCutDmn::get_elements()[k_ind], KClusterType::get_super_basis_vectors());

        for (int K_ind = 0; K_ind < KClusterDmn::dmn_size(); ++K_ind) {
          length_and_distance[K_ind].second = K_ind;

          length_and_distance[K_ind].first = domains::cluster_operations::minimal_distance(
              k_vec, KClusterDmn::get_elements()[K_ind], KClusterType::get_super_basis_vectors());
        }

        std::sort(length_and_distance.begin(), length_and_distance.end());

        int result_ind = length_and_distance[0].second;

        for (int nu_ind = 0; nu_ind < 2 * BDmn::dmn_size(); ++nu_ind) {
          Sigma_band_structure(nu_ind, k_ind) =
              Sigma(nu_ind, nu_ind, result_ind, WDmn::dmn_size() / 2);
          Sigma_cluster_band_structure(nu_ind, k_ind) =
              Sigma_cluster(nu_ind, nu_ind, result_ind, WDmn::dmn_size() / 2);
        }
      }
    }
  }

  Sigma_lattice_band_structure.reset();
  if (parameters_.do_dca_plus()) {
    func::function<std::complex<Real>, func::dmn_variadic<NuDmn, KHostDmn>> S_k_dmn("S_k_dmn_s");

    for (int b_ind = 0; b_ind < BDmn::dmn_size(); ++b_ind)
      for (int s_ind = 0; s_ind < SDmn::dmn_size(); ++s_ind)
        for (int k_ind = 0; k_ind < KHostDmn::dmn_size(); ++k_ind)
          S_k_dmn(b_ind, s_ind, k_ind) =
              Sigma_lattice(b_ind, s_ind, b_ind, s_ind, k_ind, WDmn::dmn_size() / 2);

    domains::hspline_interpolation<KHostDmn, KCutDmn>::execute(
        S_k_dmn, Sigma_lattice_band_structure, -1. / 2.);
  }

  Sigma_band_structure_interpolated.reset();

  func::function<std::complex<Real>, func::dmn_variadic<NuDmn, KHostDmn>> S_k_dmn("S_k_dmn_s");

  for (int b_ind = 0; b_ind < BDmn::dmn_size(); ++b_ind)
    for (int s_ind = 0; s_ind < SDmn::dmn_size(); ++s_ind)
      for (int k_ind = 0; k_ind < KHostDmn::dmn_size(); ++k_ind)
        S_k_dmn(b_ind, s_ind, k_ind) =
            Sigma_lattice_interpolated(b_ind, s_ind, b_ind, s_ind, k_ind, WDmn::dmn_size() / 2);

  domains::hspline_interpolation<KHostDmn, KCutDmn>::execute(
      S_k_dmn, Sigma_band_structure_interpolated, -1. / 2.);

  Sigma_band_structure_coarsegrained.reset();
  if (parameters_.do_dca_plus()) {
    func::function<std::complex<Real>, func::dmn_variadic<NuDmn, KHostDmn>> S_k_dmn("S_k_dmn_s");

    for (int b_ind = 0; b_ind < BDmn::dmn_size(); ++b_ind)
      for (int s_ind = 0; s_ind < SDmn::dmn_size(); ++s_ind)
        for (int k_ind = 0; k_ind < KHostDmn::dmn_size(); ++k_ind)
          S_k_dmn(b_ind, s_ind, k_ind) =
              Sigma_lattice_coarsegrained(b_ind, s_ind, b_ind, s_ind, k_ind, WDmn::dmn_size() / 2);

    domains::hspline_interpolation<KHostDmn, KCutDmn>::execute(
        S_k_dmn, Sigma_band_structure_coarsegrained, -1. / 2.);
  }
}

template <class Parameters, DistType DT>
void DcaData<Parameters, DT>::print_Sigma_QMC_versus_Sigma_cg() {
  if (concurrency_.id() == concurrency_.first() /*and parameters_.do_dca_plus()*/) {
    if (DIMENSION == 2) {
      std::cout << "\n\n";
      std::cout << "        K-vectors             || Re[Sigma_QMC]   Im[Sigma_QMC]   Re[Sigma_cg]  "
                   "  Im[Sigma_cg] \n";
      std::cout << "-------------------------------------------------------------------------------"
                   "---------------\n";
    }

    if (DIMENSION == 3) {
      std::cout << "\n\n";
      std::cout << "                K-vectors                       || Re[Sigma_QMC]   "
                   "Im[Sigma_QMC]   Re[Sigma_cg]    Im[Sigma_cg] \n";
      std::cout << "-------------------------------------------------------------------------------"
                   "---------------------------------\n";
    }

    for (int k_ind = 0; k_ind < KClusterDmn::dmn_size(); ++k_ind) {
      math::util::print(KClusterDmn::get_elements()[k_ind]);
      std::cout << real(Sigma(0, 0, k_ind, WDmn::dmn_size() / 2)) << "\t"
                << imag(Sigma(0, 0, k_ind, WDmn::dmn_size() / 2)) << "\t";
      std::cout << real(Sigma_cluster(0, 0, k_ind, WDmn::dmn_size() / 2)) << "\t"
                << imag(Sigma_cluster(0, 0, k_ind, WDmn::dmn_size() / 2)) << "\n";
    }
    std::cout << "\n\n";
  }
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_DATA_DCA_DATA_HPP
