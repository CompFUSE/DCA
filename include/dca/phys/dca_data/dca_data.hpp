// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
#include "dca/phys/domains/quantum/brillouin_zone_cut_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/four_point_type.hpp"
#include "dca/phys/parameters/cluster_domain_aliases.hpp"
#include "dca/util/timer.hpp"


namespace dca {
namespace phys {
// dca::phys::

template <class parameters_type>
class DcaData {
public:
  constexpr static int DIMENSION = parameters_type::lattice_type::DIMENSION;

  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using CDA = ClusterDomainAliases<DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterType = typename CDA::KClusterType;
  using KClusterDmn = typename CDA::KClusterDmn;
  using r_HOST = typename CDA::RSpHostDmn;
  using k_HOST = typename CDA::KSpHostDmn;

  using t = func::dmn_0<domains::time_domain>;
  using w = func::dmn_0<domains::frequency_domain>;
  using w_VERTEX = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index
  using nu_nu = func::dmn_variadic<nu, nu>;

  using k_domain_cut_dmn_type = func::dmn_0<domains::brillouin_zone_cut_domain<101>>;
  using nu_k_cut = func::dmn_variadic<nu, k_domain_cut_dmn_type>;

  using nu_nu_k_DCA_w = func::dmn_variadic<nu, nu, KClusterDmn, w>;


  DcaData(parameters_type& parameters_ref);

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
  parameters_type& parameters;
  concurrency_type& concurrency;

public:
  func::function<int, nu_nu> H_symmetry;
  func::function<double, func::dmn_variadic<nu, nu, RClusterDmn>> H_interactions;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn>> H_DCA;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_HOST>> H_HOST;

  func::function<double, nu_k_cut> band_structure;

  func::function<std::complex<double>, nu_k_cut> Sigma_band_structure;

  func::function<std::complex<double>, nu_k_cut> Sigma_cluster_band_structure;
  func::function<std::complex<double>, nu_k_cut> Sigma_lattice_band_structure;

  func::function<std::complex<double>, nu_k_cut> Sigma_band_structure_interpolated;
  func::function<std::complex<double>, nu_k_cut> Sigma_band_structure_coarsegrained;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_HOST>> G_k;  //("Greens-k-lattice");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_HOST>> S_k;  //("Sigma-k-lattice");
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_HOST>> S_r;  //("Sigma-r-lattice");

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> Sigma;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> Sigma_stddev;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> Sigma_cluster;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_HOST, w>> Sigma_lattice;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_HOST, w>> Sigma_lattice_interpolated;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, k_HOST, w>> Sigma_lattice_coarsegrained;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> G_k_w;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> G_k_w_stddev;
  func::function<double, func::dmn_variadic<nu, nu, KClusterDmn, t>> G_k_t;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>> G_r_w;
  func::function<double, func::dmn_variadic<nu, nu, RClusterDmn, t>> G_r_t;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> G0_k_w;
  func::function<double, func::dmn_variadic<nu, nu, KClusterDmn, t>> G0_k_t;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>> G0_r_w;
  func::function<double, func::dmn_variadic<nu, nu, RClusterDmn, t>> G0_r_t;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn, w>> G0_k_w_cluster_excluded;
  func::function<double, func::dmn_variadic<nu, nu, KClusterDmn, t>> G0_k_t_cluster_excluded;
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn, w>> G0_r_w_cluster_excluded;
  func::function<double, func::dmn_variadic<nu, nu, RClusterDmn, t>> G0_r_t_cluster_excluded;

  func::function<std::complex<double>, func::dmn_variadic<b, b, b, b, KClusterDmn, KClusterDmn, w_VERTEX, w_VERTEX>>
      G4_k_k_w_w;
  func::function<std::complex<double>, func::dmn_variadic<b, b, b, b, KClusterDmn, KClusterDmn, w_VERTEX, w_VERTEX>>
      G4_k_k_w_w_stddev;

  func::function<double, func::dmn_variadic<nu, nu, RClusterDmn, t>> K_r_t;

  func::function<double, nu> orbital_occupancy;
};

template <class parameters_type>
DcaData<parameters_type>::DcaData(parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

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
      Sigma_stddev("Self_Energy-stddev"),

      Sigma_cluster("Self-Energy-cluster"),
      Sigma_lattice("Self-energy-lattice"),

      Sigma_lattice_interpolated("Sigma_lattice_interpolated"),
      Sigma_lattice_coarsegrained("Sigma_lattice_coarsegrained"),

      G_k_w("cluster_greens_function_G_k_w"),
      G_k_w_stddev("cluster_greens_function_G_k_w-stddev"),
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

      G4_k_k_w_w("G4_k_k_w_w"),
      G4_k_k_w_w_stddev("G4_k_k_w_w-stddev"),

      K_r_t("K_r_t"),

      orbital_occupancy("orbital_occupancy") {
  H_symmetry = -1;
}

template <class parameters_type>
void DcaData<parameters_type>::read(std::string filename) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t starts reading \n\n";

  if (concurrency.id() == concurrency.first()) {
    const std::string& output_format = parameters.get_output_format();

    if (output_format == "JSON") {
      dca::io::JSONReader reader;
      reader.open_file(filename);
      this->read(reader);
      reader.close_file();
    }

    else if (output_format == "HDF5") {
      dca::io::HDF5Reader reader;
      reader.open_file(filename);
      this->read(reader);
      reader.close_file();
    }

    else
      throw std::logic_error(__FUNCTION__);
  }

  concurrency.broadcast(parameters.get_chemical_potential());

  concurrency.broadcast_object(Sigma);

  if (parameters.get_four_point_type() != NONE) {
    concurrency.broadcast_object(G_k_w);
    concurrency.broadcast_object(G4_k_k_w_w);
  }
}

template <class parameters_type>
template <typename Reader>
void DcaData<parameters_type>::read(Reader& reader) {
  std::string four_point_type = "NONE";

  {
    reader.open_group("parameters");

    {
      reader.open_group("physics");

      reader.execute("chemical-potential", parameters.get_chemical_potential());

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

    if (four_point_type != "NONE") {
      reader.execute(G_k_w);

      reader.execute(G4_k_k_w_w);
    }

    reader.close_group();
  }
}

template <class parameters_type>
void DcaData<parameters_type>::write(std::string file_name) {
  std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

  const std::string& output_format = parameters.get_output_format();

  if (output_format == "JSON") {
    dca::io::JSONWriter writer;
    writer.open_file(file_name);

    parameters.write(writer);
    this->write(writer);

    writer.close_file();
  }

  else if (output_format == "HDF5") {
    dca::io::HDF5Writer writer;
    writer.open_file(file_name);

    parameters.write(writer);
    this->write(writer);

    writer.close_file();
  }

  else
    throw std::logic_error(__FUNCTION__);
}

template <class parameters_type>
template <typename Writer>
void DcaData<parameters_type>::write(Writer& writer) {
  writer.open_group("functions");

  writer.execute(band_structure);

  if (parameters.do_dca_plus()) {
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
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn>> S_r_DCA("Sigma-r-DCA");

    func::function<std::complex<double>, func::dmn_variadic<nu, nu, KClusterDmn>> S_k_DCA("Sigma-k-DCA");
    std::memcpy(&S_k_DCA(0), &Sigma(0, 0, 0, w::dmn_size() / 2),
                sizeof(std::complex<double>) * std::pow(2 * b::dmn_size(), 2.) * KClusterDmn::dmn_size());
    math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(S_k_DCA, S_r_DCA);

    writer.execute(S_r_DCA);
  }

  writer.execute(Sigma);
  writer.execute(Sigma_stddev);

  if (parameters.dump_lattice_self_energy()) {
    writer.execute(Sigma_lattice);
  }

  if (parameters.dump_cluster_Greens_functions()) {
    writer.execute(G_k_w);
    writer.execute(G_k_w_stddev);
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

  if (parameters.get_four_point_type() != NONE) {
    if (!(parameters.dump_cluster_Greens_functions())) {
      writer.execute(G_k_w);
      writer.execute(G_k_w_stddev);
    }

    writer.execute(G4_k_k_w_w);
    writer.execute(G4_k_k_w_w_stddev);
  }

  writer.close_group();
}

template <class parameters_type>
void DcaData<parameters_type>::initialize() {
  initialize_H_0_and_H_i();
  initialize_G0();
}

template <class parameters_type>
void DcaData<parameters_type>::initialize_H_0_and_H_i() {
  util::Timer("H_0 and H_int initialization", concurrency.id() == concurrency.first());

  parameters_type::model_type::initialize_H_0(parameters, H_DCA);
  parameters_type::model_type::initialize_H_0(parameters, H_HOST);

  parameters_type::model_type::initialize_H_interaction(H_interactions, parameters);
  parameters_type::model_type::initialize_H_symmetries(H_symmetry);

  compute_band_structure::execute(parameters, band_structure);
}

template <class parameters_type>
void DcaData<parameters_type>::initialize_G0() {
  profiler_type prof(__FUNCTION__, "DcaData", __LINE__);

  util::Timer("G_0 initialization", concurrency.id() == concurrency.first());

  // Compute G0_k_w.
  compute_G0_k_w(H_DCA, parameters.get_chemical_potential(), concurrency, G0_k_w);
  symmetrize::execute(G0_k_w, H_symmetry, true);

  // Compute G0_k_t.
  compute_G0_k_t(H_DCA, parameters.get_chemical_potential(), parameters.get_beta(), G0_k_t);
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

template <class parameters_type>
void DcaData<parameters_type>::initialize_Sigma() {
  profiler_type prof("initialize-Sigma", "input", __LINE__);

  if (parameters.get_initial_self_energy() != "zero")
    this->read(parameters.get_initial_self_energy());
}

template <class parameters_type>
void DcaData<parameters_type>::compute_single_particle_properties() {
  {
    std::memcpy(&S_k(0), &Sigma_lattice(0, 0, 0, w::dmn_size() / 2),
                sizeof(std::complex<double>) * std::pow(2 * b::dmn_size(), 2.) * k_HOST::dmn_size());

    math::transform::FunctionTransform<k_HOST, r_HOST>::execute(S_k, S_r);
  }

  {
    int w_ind = w::dmn_size() / 2;

    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> I_k("I_matrix", nu::dmn_size());
    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_inv("G_inv", nu::dmn_size());

    // Allocate the work space for inverse only once.
    dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
    dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> work;

    std::complex<double> i_wm_plus_mu;

    i_wm_plus_mu.real(parameters.get_chemical_potential());
    i_wm_plus_mu.imag(w::get_elements()[w_ind]);

    for (int i = 0; i < nu::dmn_size(); i++)
      I_k(i, i) = i_wm_plus_mu;

    for (int k_ind = 0; k_ind < k_HOST::dmn_size(); k_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_inv(i, j) = I_k(i, j) - H_HOST(i, j, k_ind) - Sigma_lattice(i, j, k_ind, w_ind);

      dca::linalg::matrixop::inverse(G_inv, ipiv, work);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_k(i, j, k_ind) = G_inv(i, j);
    }
  }
}

template <class parameters_type>
void DcaData<parameters_type>::compute_Sigma_bands() {
  {
    Sigma_band_structure.reset();
    Sigma_cluster_band_structure.reset();

    std::vector<std::pair<double, int>> length_and_distance(KClusterDmn::dmn_size(),
                                                            std::pair<double, int>(0, -1));

    for (int k_ind = 0; k_ind < k_domain_cut_dmn_type::dmn_size(); ++k_ind) {
      std::vector<double> k_vec = domains::cluster_operations::translate_inside_cluster(
          k_domain_cut_dmn_type::get_elements()[k_ind],
          KClusterType::get_super_basis_vectors());

      for (int K_ind = 0; K_ind < KClusterDmn::dmn_size(); ++K_ind) {
        length_and_distance[K_ind].second = K_ind;

        length_and_distance[K_ind].first = domains::cluster_operations::minimal_distance(
            k_vec, KClusterDmn::get_elements()[K_ind], KClusterType::get_super_basis_vectors());
      }

      std::sort(length_and_distance.begin(), length_and_distance.end());

      int result_ind = length_and_distance[0].second;

      for (int nu_ind = 0; nu_ind < 2 * b::dmn_size(); ++nu_ind) {
        Sigma_band_structure(nu_ind, k_ind) = Sigma(nu_ind, nu_ind, result_ind, w::dmn_size() / 2);
        Sigma_cluster_band_structure(nu_ind, k_ind) =
            Sigma_cluster(nu_ind, nu_ind, result_ind, w::dmn_size() / 2);
      }
    }
  }

  Sigma_lattice_band_structure.reset();
  if (parameters.do_dca_plus()) {
    func::function<std::complex<double>, func::dmn_variadic<nu, k_HOST>> S_k_dmn("S_k_dmn_s");

    for (int b_ind = 0; b_ind < b::dmn_size(); ++b_ind)
      for (int s_ind = 0; s_ind < s::dmn_size(); ++s_ind)
        for (int k_ind = 0; k_ind < k_HOST::dmn_size(); ++k_ind)
          S_k_dmn(b_ind, s_ind, k_ind) =
              Sigma_lattice(b_ind, s_ind, b_ind, s_ind, k_ind, w::dmn_size() / 2);

    domains::hspline_interpolation<k_HOST, k_domain_cut_dmn_type>::execute(
        S_k_dmn, Sigma_lattice_band_structure, -1. / 2.);
  }

  Sigma_band_structure_interpolated.reset();
  if (true) {
    func::function<std::complex<double>, func::dmn_variadic<nu, k_HOST>> S_k_dmn("S_k_dmn_s");

    for (int b_ind = 0; b_ind < b::dmn_size(); ++b_ind)
      for (int s_ind = 0; s_ind < s::dmn_size(); ++s_ind)
        for (int k_ind = 0; k_ind < k_HOST::dmn_size(); ++k_ind)
          S_k_dmn(b_ind, s_ind, k_ind) =
              Sigma_lattice_interpolated(b_ind, s_ind, b_ind, s_ind, k_ind, w::dmn_size() / 2);

    domains::hspline_interpolation<k_HOST, k_domain_cut_dmn_type>::execute(
        S_k_dmn, Sigma_band_structure_interpolated, -1. / 2.);
  }

  Sigma_band_structure_coarsegrained.reset();
  if (parameters.do_dca_plus()) {
    func::function<std::complex<double>, func::dmn_variadic<nu, k_HOST>> S_k_dmn("S_k_dmn_s");

    for (int b_ind = 0; b_ind < b::dmn_size(); ++b_ind)
      for (int s_ind = 0; s_ind < s::dmn_size(); ++s_ind)
        for (int k_ind = 0; k_ind < k_HOST::dmn_size(); ++k_ind)
          S_k_dmn(b_ind, s_ind, k_ind) =
              Sigma_lattice_coarsegrained(b_ind, s_ind, b_ind, s_ind, k_ind, w::dmn_size() / 2);

    domains::hspline_interpolation<k_HOST, k_domain_cut_dmn_type>::execute(
        S_k_dmn, Sigma_band_structure_coarsegrained, -1. / 2.);
  }
}

template <class parameters_type>
void DcaData<parameters_type>::print_Sigma_QMC_versus_Sigma_cg() {
  if (concurrency.id() == concurrency.first() /*and parameters.do_dca_plus()*/) {
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
      std::cout << real(Sigma(0, 0, k_ind, w::dmn_size() / 2)) << "\t"
                << imag(Sigma(0, 0, k_ind, w::dmn_size() / 2)) << "\t";
      std::cout << real(Sigma_cluster(0, 0, k_ind, w::dmn_size() / 2)) << "\t"
                << imag(Sigma_cluster(0, 0, k_ind, w::dmn_size() / 2)) << "\n";
    }
    std::cout << "\n\n";
  }
}

}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_DATA_DCA_DATA_HPP
