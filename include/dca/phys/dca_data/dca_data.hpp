// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class contains all functions needed for the MOMS DCA calculation.
//
// TODO: Interpolate Sigma if \beta from file != \beta ?

#ifndef DCA_PHYS_DCA_DATA_DCA_DATA_HPP
#define DCA_PHYS_DCA_DATA_DCA_DATA_HPP

#include <algorithm>
#include <complex>
#include <cstring>
#include <iostream>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/io/json/json_writer.hpp"
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
#include "dca/phys/four_point_type.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/models/traits.hpp"
#include "dca/util/timer.hpp"

namespace dca {
namespace phys {
// dca::phys::

template <class Parameters>
class DcaData {
public:
  using profiler_type = typename Parameters::profiler_type;

  using Concurrency = typename Parameters::concurrency_type;
  using Lattice = typename Parameters::lattice_type;
  constexpr static int DIMENSION = Lattice::DIMENSION;
  using TpAccumulatorScalar = typename Parameters::MC_measurement_scalar_type;

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
      func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>>;
  using TpGreensFunction =
      func::function<std::complex<TpAccumulatorScalar>,
                     func::dmn_variadic<BDmn, BDmn, BDmn, BDmn, KClusterDmn, KClusterDmn,
                                        KExchangeDmn, WVertexDmn, WVertexDmn, WExchangeDmn>>;
  // The following typedef is for testing purposes
  using ReducedTpGreensFunction = func::function<
      std::complex<TpAccumulatorScalar>,
      func::dmn_variadic<BDmn, BDmn, BDmn, BDmn, KClusterDmn, KClusterDmn, WVertexDmn, WVertexDmn>>;

  DcaData(Parameters& parameters_ref);

  void read(std::string filename);
  template <typename Reader>
  void read(Reader& reader);

  void write(std::string filename);
  template <typename Writer>
  void write(Writer& reader);

  void initialize();
  void initialize_H_0_and_H_i();
  void initialize_G0();
  void initialize_Sigma();

  void compute_single_particle_properties();
  void compute_Sigma_bands();

  void print_Sigma_QMC_versus_Sigma_cg();

private:
  Parameters& parameters_;
  Concurrency& concurrency_;

public:
  func::function<int, NuNuDmn> H_symmetry;
  func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn>> H_interactions;

  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn>> H_DCA;
  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn>> H_HOST;

  func::function<double, NuKCutDmn> band_structure;

  func::function<std::complex<double>, NuKCutDmn> Sigma_band_structure;

  func::function<std::complex<double>, NuKCutDmn> Sigma_cluster_band_structure;
  func::function<std::complex<double>, NuKCutDmn> Sigma_lattice_band_structure;

  func::function<std::complex<double>, NuKCutDmn> Sigma_band_structure_interpolated;
  func::function<std::complex<double>, NuKCutDmn> Sigma_band_structure_coarsegrained;

  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn>> G_k;  //("Greens-k-lattice");
  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn>> S_k;  //("Sigma-k-lattice");
  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RHostDmn>> S_r;  //("Sigma-r-lattice");

  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>> Sigma;

  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>> Sigma_cluster;
  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn, WDmn>> Sigma_lattice;
  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn, WDmn>>
      Sigma_lattice_interpolated;
  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KHostDmn, WDmn>>
      Sigma_lattice_coarsegrained;

  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>> G_k_w;
  func::function<double, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, TDmn>> G_k_t;
  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>> G_r_w;
  func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>> G_r_t;

  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>> G0_k_w;
  func::function<double, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, TDmn>> G0_k_t;
  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>> G0_r_w;
  func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>> G0_r_t;

  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>>
      G0_k_w_cluster_excluded;
  func::function<double, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, TDmn>> G0_k_t_cluster_excluded;
  func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, WDmn>>
      G0_r_w_cluster_excluded;
  func::function<double, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn, TDmn>> G0_r_t_cluster_excluded;

  func::function<double, NuDmn> orbital_occupancy;

public:  // Optional members getters.
  auto& get_G_k_w_error() {
    if (not G_k_w_err_)
      G_k_w_err_.reset(new SpGreensFunction("G_k_w-error"));
    return *G_k_w_err_;
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
    if (not G4_)
      G4_.reset(new TpGreensFunction("G4"));
    return *G4_;
  }
  auto& get_G4_error() {
    if (not G4_err_)
      G4_err_.reset(new TpGreensFunction("G4-error"));
    return *G4_err_;
  }
  auto& get_G4_stdv() {
    if (not G4_err_)
      G4_err_.reset(new TpGreensFunction("G4-stddev"));
    return *G4_err_;
  }
  auto& get_non_density_interactions() {
    if (not non_density_interactions_)
      non_density_interactions_.reset(
          new func::function<double, func::dmn_variadic<NuDmn, NuDmn, NuDmn, NuDmn, RClusterDmn>>(
              "non_density_interaction"));
    return *non_density_interactions_;
  }
  bool has_non_density_interactions() const {
    return (bool)non_density_interactions_;
  }

private:  // Optional members.
  std::unique_ptr<SpGreensFunction> G_k_w_err_;
  std::unique_ptr<SpGreensFunction> Sigma_err_;
  std::unique_ptr<TpGreensFunction> G4_;
  std::unique_ptr<TpGreensFunction> G4_err_;
  std::unique_ptr<func::function<double, func::dmn_variadic<NuDmn, NuDmn, NuDmn, NuDmn, RClusterDmn>>>
      non_density_interactions_;
};

template <class Parameters>
DcaData<Parameters>::DcaData(Parameters& parameters_ref)
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
}

template <class Parameters>
void DcaData<Parameters>::read(std::string filename) {
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\n\t starts reading \n\n";

  if (concurrency_.id() == concurrency_.first()) {
    const std::string& output_format = parameters_.get_output_format();

    if (output_format == static_cast<const std::string>("JSON")) {
      dca::io::JSONReader reader;
      reader.open_file(filename);
      this->read(reader);
      reader.close_file();
    }

    else if (output_format == static_cast<const std::string>("HDF5")) {
      dca::io::HDF5Reader reader;
      reader.open_file(filename);
      this->read(reader);
      reader.close_file();
    }

    else
      throw std::logic_error(__FUNCTION__);
  }

  concurrency_.broadcast(parameters_.get_chemical_potential());

  concurrency_.broadcast_object(Sigma);

  if (parameters_.get_four_point_type() != NONE) {
    concurrency_.broadcast_object(G_k_w);
    concurrency_.broadcast_object(get_G4());
  }
}

template <class Parameters>
template <typename Reader>
void DcaData<Parameters>::read(Reader& reader) {
  std::string four_point_type = "NONE";

  {
    reader.open_group("parameters");

    {
      reader.open_group("physics");

      reader.execute("chemical-potential", parameters_.get_chemical_potential());

      reader.close_group();
    }

    {
      reader.open_group("four-point");

      reader.execute("type", four_point_type);

      reader.close_group();
    }

    reader.close_group();
  }

  {
    reader.open_group("functions");

    reader.execute(Sigma);

    if (four_point_type != static_cast<const std::string>("NONE")) {
      reader.execute(G_k_w);

      reader.execute(get_G4());
    }

    reader.close_group();
  }
}

template <class Parameters>
void DcaData<Parameters>::write(std::string file_name) {
  std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

  const std::string& output_format = parameters_.get_output_format();

  if (output_format == static_cast<const std::string>("JSON")) {
    dca::io::JSONWriter writer;
    writer.open_file(file_name);

    parameters_.write(writer);
    this->write(writer);

    writer.close_file();
  }

  else if (output_format == static_cast<const std::string>("HDF5")) {
    dca::io::HDF5Writer writer;
    writer.open_file(file_name);

    parameters_.write(writer);
    this->write(writer);

    writer.close_file();
  }

  else
    throw std::logic_error(__FUNCTION__);
}

template <class Parameters>
template <typename Writer>
void DcaData<Parameters>::write(Writer& writer) {
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
    func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, RClusterDmn>> S_r_DCA(
        "Sigma-r-DCA");

    func::function<std::complex<double>, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn>> S_k_DCA(
        "Sigma-k-DCA");
    std::memcpy(&S_k_DCA(0), &Sigma(0, 0, 0, WDmn::dmn_size() / 2),
                sizeof(std::complex<double>) * std::pow(2 * BDmn::dmn_size(), 2.) *
                    KClusterDmn::dmn_size());
    math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(S_k_DCA, S_r_DCA);

    writer.execute(S_r_DCA);
  }

  writer.execute(Sigma);
  writer.execute(Sigma_err_);

  if (parameters_.dump_lattice_self_energy()) {
    writer.execute(Sigma_lattice);
  }

  if (parameters_.dump_cluster_Greens_functions()) {
    writer.execute(G_k_w);
    writer.execute(G_k_w_err_);
    writer.execute(G_r_w);
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

  if (parameters_.get_four_point_type() != NONE) {
    if (!(parameters_.dump_cluster_Greens_functions())) {
      writer.execute(G_k_w);
      writer.execute(G_k_w_err_);
    }

    writer.execute(G4_);
    writer.execute(G4_err_);
  }

  writer.close_group();
}

template <class Parameters>
void DcaData<Parameters>::initialize() {
  initialize_H_0_and_H_i();
  initialize_G0();
}

template <class Parameters>
void DcaData<Parameters>::initialize_H_0_and_H_i() {
  util::Timer("H_0 and H_int initialization", concurrency_.id() == concurrency_.first());

  Parameters::model_type::initialize_H_0(parameters_, H_DCA);
  Parameters::model_type::initialize_H_0(parameters_, H_HOST);

  Parameters::model_type::initialize_H_interaction(H_interactions, parameters_);

  if (models::has_non_density_interaction<Lattice>::value) {
    models::initializeNonDensityInteraction<Lattice>(get_non_density_interactions(), parameters_);
  }

  Parameters::model_type::initialize_H_symmetries(H_symmetry);

  compute_band_structure::execute(parameters_, band_structure);
}

template <class Parameters>
void DcaData<Parameters>::initialize_G0() {
  profiler_type prof(__FUNCTION__, "DcaData", __LINE__);

  util::Timer("G_0 initialization", concurrency_.id() == concurrency_.first());

  // Compute G0_k_w.
  compute_G0_k_w(H_DCA, parameters_.get_chemical_potential(), concurrency_, G0_k_w);
  symmetrize::execute(G0_k_w, H_symmetry, true);

  // Compute G0_k_t.
  compute_G0_k_t(H_DCA, parameters_.get_chemical_potential(), parameters_.get_beta(), G0_k_t);
  symmetrize::execute(G0_k_t, H_symmetry, true);

  // Compute G0_r_w.
  math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(G0_k_w, G0_r_w);
  symmetrize::execute(G0_r_w, H_symmetry, true);

  // Compute G0_r_t.
  math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(G0_k_t, G0_r_t);
  symmetrize::execute(G0_r_t, H_symmetry, true);

  // Initialize the cluster excluded Green's functions with the corresponding free Green's
  // functions.
  G0_k_w_cluster_excluded = G0_k_w;
  G0_k_t_cluster_excluded = G0_k_t;
  G0_r_w_cluster_excluded = G0_r_w;
  G0_r_t_cluster_excluded = G0_r_t;
}

template <class Parameters>
void DcaData<Parameters>::initialize_Sigma() {
  profiler_type prof("initialize-Sigma", "input", __LINE__);

  if (parameters_.get_initial_self_energy() != "zero")
    this->read(parameters_.get_initial_self_energy());
}

template <class Parameters>
void DcaData<Parameters>::compute_single_particle_properties() {
  {
    std::memcpy(
        &S_k(0), &Sigma_lattice(0, 0, 0, WDmn::dmn_size() / 2),
        sizeof(std::complex<double>) * std::pow(2 * BDmn::dmn_size(), 2.) * KHostDmn::dmn_size());

    math::transform::FunctionTransform<KHostDmn, RHostDmn>::execute(S_k, S_r);
  }

  {
    int w_ind = WDmn::dmn_size() / 2;

    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> I_k("I_matrix", NuDmn::dmn_size());
    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_inv("G_inv", NuDmn::dmn_size());

    // Allocate the work space for inverse only once.
    dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
    dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> work;

    std::complex<double> i_wm_plus_mu;

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

template <class Parameters>
void DcaData<Parameters>::compute_Sigma_bands() {
  {
    Sigma_band_structure.reset();
    Sigma_cluster_band_structure.reset();

    std::vector<std::pair<double, int>> length_and_distance(KClusterDmn::dmn_size(),
                                                            std::pair<double, int>(0, -1));

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
        Sigma_band_structure(nu_ind, k_ind) = Sigma(nu_ind, nu_ind, result_ind, WDmn::dmn_size() / 2);
        Sigma_cluster_band_structure(nu_ind, k_ind) =
            Sigma_cluster(nu_ind, nu_ind, result_ind, WDmn::dmn_size() / 2);
      }
    }
  }

  Sigma_lattice_band_structure.reset();
  if (parameters_.do_dca_plus()) {
    func::function<std::complex<double>, func::dmn_variadic<NuDmn, KHostDmn>> S_k_dmn("S_k_dmn_s");

    for (int b_ind = 0; b_ind < BDmn::dmn_size(); ++b_ind)
      for (int s_ind = 0; s_ind < SDmn::dmn_size(); ++s_ind)
        for (int k_ind = 0; k_ind < KHostDmn::dmn_size(); ++k_ind)
          S_k_dmn(b_ind, s_ind, k_ind) =
              Sigma_lattice(b_ind, s_ind, b_ind, s_ind, k_ind, WDmn::dmn_size() / 2);

    domains::hspline_interpolation<KHostDmn, KCutDmn>::execute(
        S_k_dmn, Sigma_lattice_band_structure, -1. / 2.);
  }

  Sigma_band_structure_interpolated.reset();

  func::function<std::complex<double>, func::dmn_variadic<NuDmn, KHostDmn>> S_k_dmn("S_k_dmn_s");

  for (int b_ind = 0; b_ind < BDmn::dmn_size(); ++b_ind)
    for (int s_ind = 0; s_ind < SDmn::dmn_size(); ++s_ind)
      for (int k_ind = 0; k_ind < KHostDmn::dmn_size(); ++k_ind)
        S_k_dmn(b_ind, s_ind, k_ind) =
            Sigma_lattice_interpolated(b_ind, s_ind, b_ind, s_ind, k_ind, WDmn::dmn_size() / 2);

  domains::hspline_interpolation<KHostDmn, KCutDmn>::execute(
      S_k_dmn, Sigma_band_structure_interpolated, -1. / 2.);

  Sigma_band_structure_coarsegrained.reset();
  if (parameters_.do_dca_plus()) {
    func::function<std::complex<double>, func::dmn_variadic<NuDmn, KHostDmn>> S_k_dmn("S_k_dmn_s");

    for (int b_ind = 0; b_ind < BDmn::dmn_size(); ++b_ind)
      for (int s_ind = 0; s_ind < SDmn::dmn_size(); ++s_ind)
        for (int k_ind = 0; k_ind < KHostDmn::dmn_size(); ++k_ind)
          S_k_dmn(b_ind, s_ind, k_ind) =
              Sigma_lattice_coarsegrained(b_ind, s_ind, b_ind, s_ind, k_ind, WDmn::dmn_size() / 2);

    domains::hspline_interpolation<KHostDmn, KCutDmn>::execute(
        S_k_dmn, Sigma_band_structure_coarsegrained, -1. / 2.);
  }
}

template <class Parameters>
void DcaData<Parameters>::print_Sigma_QMC_versus_Sigma_cg() {
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

}  //  phys
}  // dca

#endif  // DCA_PHYS_DCA_DATA_DCA_DATA_HPP
