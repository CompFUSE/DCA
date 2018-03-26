// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
// See LICENSE.txt for terms of usage./
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DATA_CT_INT_DATA_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DATA_CT_INT_DATA_HPP

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
#include "dca/function/function_utils.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/io/json/json_writer.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/dca_algorithms/compute_band_structure.hpp"
#include "dca/phys/dca_algorithms/compute_free_greens_function.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_sp.hpp"
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
#include "dca/phys/models/traits.hpp"
#include "dca/util/timer.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <class parameters_type>
class DcaDataCtInt {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using Lattice = typename parameters_type::lattice_type;

  using Tdmn = func::dmn_0<domains::time_domain>;
  using Wdmn = func::dmn_0<domains::frequency_domain>;

  using Bdmn = func::dmn_0<domains::electron_band_domain>;
  using Sdmn = func::dmn_0<domains::electron_spin_domain>;
  using Nu = func::dmn_variadic<Bdmn, Sdmn>;  // orbital-spin index
  using NuNu = func::dmn_variadic<Nu, Nu>;

  constexpr static int DIMENSION = parameters_type::lattice_type::DIMENSION;
  using Rdmn = func::dmn_0<domains::cluster_domain<double, DIMENSION, domains::CLUSTER,
                                                   domains::REAL_SPACE, domains::BRILLOUIN_ZONE>>;
  using Kdmn = func::dmn_0<domains::cluster_domain<double, DIMENSION, domains::CLUSTER,
                                                   domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  using k_HOST =
      func::dmn_0<domains::cluster_domain<double, DIMENSION, domains::LATTICE_SP,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;

  using k_domain_cut_dmn_type = func::dmn_0<domains::brillouin_zone_cut_domain<101>>;
  using nu_k_cut = func::dmn_variadic<Nu, k_domain_cut_dmn_type>;
  using nu_nu_k_DCA_w = func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>;

  DcaDataCtInt(parameters_type& parameters_ref);

  void read(std::string filename);
  template <typename Reader>
  void read(Reader& reader);

  void write(const std::string& filename);
  void write(const char* filename) {
    write(std::string(filename));
  }
  template <typename Writer>
  void write(Writer& reader);

  void initialize();
  void initialize_H_0_and_H_i();
  void initialize_G0();
  void initialize_Sigma();

  // Empty method for compatibility with DcaLoop.
  void compute_single_particle_properties() {}
  void compute_Sigma_bands() {}

  auto get_H_DCA() const {
    return func::utils::complex(H_0);
  }

  void print_Sigma_QMC_versus_Sigma_cg();

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  using W2pDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using TpGreenFunction =
      func::function<std::complex<double>,
                     func::dmn_variadic<Bdmn, Bdmn, Bdmn, Bdmn, Kdmn, Kdmn, W2pDmn, W2pDmn>>;

public:
  func::function<int, NuNu> H_symmetry;

  func::function<double, func::dmn_variadic<Nu, Nu, Kdmn>> H_0;
  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, k_HOST>> H_HOST;
  func::function<double, func::dmn_variadic<Nu, Nu, Rdmn>> H_interactions;
  std::unique_ptr<func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, Rdmn>>> non_density_interactions;

  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>> Sigma;
  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>> Sigma_cluster;
  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, k_HOST, Wdmn>> Sigma_lattice;
  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, k_HOST, Wdmn>> Sigma_lattice_interpolated;
  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, k_HOST, Wdmn>> Sigma_lattice_coarsegrained;

  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>> G_k_w;
  func::function<double, func::dmn_variadic<Nu, Nu, Kdmn, Tdmn>> G_k_t;
  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Rdmn, Wdmn>> G_r_w;
  func::function<double, func::dmn_variadic<Nu, Nu, Rdmn, Tdmn>> G_r_t;

  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>> M;
  std::unique_ptr<func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>>> G_k_w_error;
  std::unique_ptr<func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>>> Sigma_error;
  std::unique_ptr<func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>>> M_error;

  std::unique_ptr<TpGreenFunction> G4_k_k_w_w;
  std::unique_ptr<TpGreenFunction> G4_k_k_w_w_error;

  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>> G0_k_w;
  func::function<double, func::dmn_variadic<Nu, Nu, Kdmn, Tdmn>> G0_k_t;
  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Rdmn, Wdmn>> G0_r_w;
  func::function<double, func::dmn_variadic<Nu, Nu, Rdmn, Tdmn>> G0_r_t;

  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>> G0_k_w_cluster_excluded;
  func::function<double, func::dmn_variadic<Nu, Nu, Kdmn, Tdmn>> G0_k_t_cluster_excluded;
  func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Rdmn, Wdmn>> G0_r_w_cluster_excluded;
  func::function<double, func::dmn_variadic<Nu, Nu, Rdmn, Tdmn>> G0_r_t_cluster_excluded;

  func::function<double, Nu> orbital_occupancy;
};

template <class parameters_type>
DcaDataCtInt<parameters_type>::DcaDataCtInt(parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      H_symmetry("H_symmetry"),

      Sigma("Self_Energy"),

      Sigma_cluster("Self-Energy-cluster"),
      Sigma_lattice("Self-energy-lattice"),

      Sigma_lattice_interpolated("Sigma_lattice_interpolated"),
      Sigma_lattice_coarsegrained("Sigma_lattice_coarsegrained"),

      G_k_w("cluster_greens_function_G_k_w"),
      G_k_t("cluster_greens_function_G_k_t"),
      G_r_w("cluster_greens_function_G_r_w"),
      G_r_t("cluster_greens_function_G_r_t"),

      M("M"),

      G0_k_w("free_cluster_greens_function_G0_k_w"),
      G0_k_t("free_cluster_greens_function_G0_k_t"),
      G0_r_w("free_cluster_greens_function_G0_r_w"),
      G0_r_t("free_cluster_greens_function_G0_r_t"),

      G0_k_w_cluster_excluded("cluster_excluded_greens_function_G0_k_w"),
      G0_k_t_cluster_excluded("cluster_excluded_greens_function_G0_k_t"),
      G0_r_w_cluster_excluded("cluster_excluded_greens_function_G0_r_w"),
      G0_r_t_cluster_excluded("cluster_excluded_greens_function_G0_r_t"),

      orbital_occupancy("orbital_occupancy") {
  if (parameters.get_four_point_type() != NONE)
    G4_k_k_w_w.reset(new TpGreenFunction("G4_k_k_w_w"));

  if (parameters.computeError()) {
    G_k_w_error.reset(new func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>>(
        "G_k_w-error"));
    Sigma_error.reset(new func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>>(
        "Sigma-error"));
    M_error.reset(new func::function<std::complex<double>, func::dmn_variadic<Nu, Nu, Kdmn, Wdmn>>(
        "M-error"));
    if (parameters.get_four_point_type() != NONE)
      G4_k_k_w_w_error.reset(new TpGreenFunction("G4_k_k_w_w-error"));
  }
  H_symmetry = -1;
}

template <class parameters_type>
void DcaDataCtInt<parameters_type>::read(std::string filename) {
  if (concurrency.id() == 0)
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

  concurrency.broadcast_object(G_k_t);
}

template <class parameters_type>
template <typename Reader>
void DcaDataCtInt<parameters_type>::read(Reader& reader) {
  std::string vertex_measurement = "NONE";

  {
    reader.open_group("parameters");

    {
      reader.open_group("physics-parameters");

      reader.execute("chemical-potential", parameters.get_chemical_potential());

      reader.close_group();
    }

    {
      reader.open_group("vertex-channel");

      reader.execute("vertex-measurement-type", vertex_measurement);

      reader.close_group();
    }

    reader.close_group();
  }

  {
    reader.open_group("functions");

    reader.execute(Sigma);

    reader.close_group();
  }
}

template <class parameters_type>
void DcaDataCtInt<parameters_type>::write(const std::string& file_name) {
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
    dca::io::HDF5Writer writer(/*verbose =*/false);
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
void DcaDataCtInt<parameters_type>::write(Writer& writer) {
  writer.open_group("functions");

  writer.execute(Sigma);

  writer.execute(G_k_w);
  writer.execute(G_r_w);
  writer.execute(G_r_t);
  writer.execute(G_k_t);

  writer.execute(G0_k_w);
  writer.execute(G0_r_w);
  writer.execute(G0_k_t);
  writer.execute(G0_r_t);

  writer.execute(G0_k_w_cluster_excluded);
  writer.execute(G0_r_w_cluster_excluded);
  writer.execute(G0_k_t_cluster_excluded);
  writer.execute(G0_r_t_cluster_excluded);

  writer.execute(M);
  writer.execute(G_k_w_error);
  writer.execute(Sigma_error);
  writer.execute(M_error);

  writer.execute(G4_k_k_w_w);
  writer.execute(G4_k_k_w_w_error);

  writer.close_group();
}

template <class parameters_type>
void DcaDataCtInt<parameters_type>::initialize() {
  initialize_H_0_and_H_i();
  initialize_G0();

  using Lattice = typename parameters_type::lattice_type;
  Lattice::initialize_H_0(parameters, H_HOST);
}

template <class parameters_type>
void DcaDataCtInt<parameters_type>::initialize_H_0_and_H_i() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t initialize H_0(k) and H_i " << dca::util::print_time() << "\n";

  Lattice::initialize_H_0(parameters, H_0);

  Lattice::initialize_H_interaction(H_interactions, parameters);

  if (models::has_non_density_interaction<Lattice>::value) {
    non_density_interactions.reset(new func::function<double, func::dmn_variadic<Nu, Nu, Nu, Nu, Rdmn>>(
        "non_density_interactions"));
    models::initializeNonDensityInteraction<Lattice>(*non_density_interactions, parameters);
  }

  parameters_type::model_type::initialize_H_symmetries(H_symmetry);

  if (concurrency.id() == concurrency.first())
    std::cout << "\t finished H_0(k) and H_i " << dca::util::print_time() << "\n";
}

template <class parameters_type>
void DcaDataCtInt<parameters_type>::initialize_G0() {
  profiler_type prof(__FUNCTION__, "DcaData", __LINE__);

  dca::util::Timer("G_0 initialization", concurrency.id() == concurrency.first());

  auto H_cmplx = func::utils::complex(H_0);

  // Compute G0_k_w.
  compute_G0_k_w(H_cmplx, parameters.get_chemical_potential(), concurrency, G0_k_w);
  symmetrize::execute(G0_k_w, H_symmetry, true);

  // Compute G0_k_t.
  compute_G0_k_t(H_cmplx, parameters.get_chemical_potential(), parameters.get_beta(), G0_k_t);
  symmetrize::execute(G0_k_t, H_symmetry, true);

  // Compute G0_r_w.
  math::transform::FunctionTransform<Kdmn, Rdmn>::execute(G0_k_w, G0_r_w);
  symmetrize::execute(G0_r_w, H_symmetry, true);

  // Compute G0_r_t.
  math::transform::FunctionTransform<Kdmn, Rdmn>::execute(G0_k_t, G0_r_t);
  symmetrize::execute(G0_r_t, H_symmetry, true);

  // Copy to cluster excluded functions.
  G0_k_w_cluster_excluded = G0_k_w;
  G0_k_t_cluster_excluded = G0_k_t;
  G0_r_w_cluster_excluded = G0_r_w;
  G0_r_t_cluster_excluded = G0_r_t;
}

template <class parameters_type>
void DcaDataCtInt<parameters_type>::initialize_Sigma() {
  profiler_type prof("initialize-Sigma", "input", __LINE__);

  if (parameters.get_initial_self_energy() != "zero")
    this->read(parameters.get_initial_self_energy());
}

template <class parameters_type>
void DcaDataCtInt<parameters_type>::print_Sigma_QMC_versus_Sigma_cg() {
  if (concurrency.id() == 0 /*and parameters.do_DCA_plus()*/) {
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

    for (int k_ind = 0; k_ind < Kdmn::dmn_size(); ++k_ind) {
      math::util::print(Kdmn::get_elements()[k_ind]);
      std::cout << real(Sigma(0, 0, k_ind, Wdmn::dmn_size() / 2)) << "\t"
                << imag(Sigma(0, 0, k_ind, Wdmn::dmn_size() / 2)) << "\t";
      std::cout << real(Sigma_cluster(0, 0, k_ind, Wdmn::dmn_size() / 2)) << "\t"
                << imag(Sigma_cluster(0, 0, k_ind, Wdmn::dmn_size() / 2)) << "\n";
    }
    std::cout << "\n\n";
  }
}

}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DATA_CT_INT_DATA_HPP
