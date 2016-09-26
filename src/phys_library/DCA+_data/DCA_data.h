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

#ifndef PHYS_LIBRARY_DCA_DATA_DCA_DATA_H
#define PHYS_LIBRARY_DCA_DATA_DCA_DATA_H

#include <algorithm>
#include <complex>
#include <cstring>
#include <iostream>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/phys/dca_algorithms/compute_band_structure.hpp"
#include "dca/util/print_time.hpp"
#include "comp_library/function_library/include_function_library.h"
#include "comp_library/linalg/linalg.hpp"
#include "comp_library/IO_library/IO.hpp"
#include "math_library/functional_transforms/function_transforms/function_transforms.hpp"
#include "math_library/geometry_library/vector_operations/vector_operations.hpp"
#include "phys_library/DCA+_step/cluster_mapping/coarsegraining_step/coarsegraining_sp.h"
#include "phys_library/DCA+_step/symmetrization/symmetrize.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/cluster/interpolation/hspline_interpolation/hspline_interpolation.hpp"
#include "phys_library/domains/cluster/interpolation/wannier_interpolation/wannier_interpolation.hpp"
#include "phys_library/domains/Quantum_domain/brillouin_zone_cut_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_compact.h"
#include "phys_library/domains/time_and_frequency/time_domain.h"
#include "phys_library/vertex_measurement_type.hpp"

using namespace dca::phys;

namespace DCA {

template <class parameters_type>
class DCA_data {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using t = dmn_0<time_domain>;
  using w = dmn_0<frequency_domain>;
  using compact_vertex_frequency_domain_type = DCA::vertex_frequency_domain<DCA::COMPACT>;
  using w_VERTEX = dmn_0<compact_vertex_frequency_domain_type>;

  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index
  using nu_nu = dmn_variadic<nu, nu>;

  using r_DCA = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, CLUSTER,
                                     REAL_SPACE, BRILLOUIN_ZONE>>;
  using DCA_k_cluster_type = cluster_domain<double, parameters_type::lattice_type::DIMENSION,
                                            CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>;
  using k_DCA = dmn_0<DCA_k_cluster_type>;
  using r_HOST = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, LATTICE_SP,
                                      REAL_SPACE, BRILLOUIN_ZONE>>;
  using k_HOST = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, LATTICE_SP,
                                      MOMENTUM_SPACE, BRILLOUIN_ZONE>>;
  using k_LDA = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, LATTICE_SP,
                                     MOMENTUM_SPACE, PARALLELLEPIPEDUM>>;

  using k_domain_cut_dmn_type = dmn_0<brillouin_zone_cut_domain<101>>;
  using nu_k_cut = dmn_variadic<nu, k_domain_cut_dmn_type>;

  using nu_nu_k_DCA_w = dmn_variadic<nu, nu, k_DCA, w>;

  const static int DIMENSION = parameters_type::lattice_type::DIMENSION;

public:
  DCA_data(parameters_type& parameters_ref);

  void read(std::string filename);

  void write(std::string filename);

  template <IO::FORMAT DATA_FORMAT>
  void read(IO::reader<DATA_FORMAT>& reader);

  template <IO::FORMAT DATA_FORMAT>
  void write(IO::writer<DATA_FORMAT>& reader);

  void initialize();

  void initialize_H_0_and_H_i();

  void initialize_G0();

  bool test_initialize_G0();

  void initialize_Sigma();

  void compute_Sigma_bands();
  void compute_single_particle_properties();

  void print_Sigma_QMC_versus_Sigma_cg();

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

public:
  FUNC_LIB::function<int, nu_nu> H_symmetry;
  FUNC_LIB::function<double, dmn_3<nu, nu, r_DCA>> H_interactions;

  FUNC_LIB::function<std::complex<double>, dmn_3<nu, nu, k_DCA>> H_DCA;
  FUNC_LIB::function<std::complex<double>, dmn_3<nu, nu, k_HOST>> H_HOST;
  FUNC_LIB::function<std::complex<double>, dmn_3<nu, nu, k_LDA>> H_LDA;

  FUNC_LIB::function<double, nu_k_cut> band_structure;

  FUNC_LIB::function<std::complex<double>, nu_k_cut> Sigma_band_structure;

  FUNC_LIB::function<std::complex<double>, nu_k_cut> Sigma_cluster_band_structure;
  FUNC_LIB::function<std::complex<double>, nu_k_cut> Sigma_lattice_band_structure;

  FUNC_LIB::function<std::complex<double>, nu_k_cut> Sigma_band_structure_interpolated;
  FUNC_LIB::function<std::complex<double>, nu_k_cut> Sigma_band_structure_coarsegrained;

  FUNC_LIB::function<std::complex<double>, dmn_3<nu, nu, k_HOST>> G_k;  //("Greens-k-lattice");
  FUNC_LIB::function<std::complex<double>, dmn_3<nu, nu, k_HOST>> S_k;  //("Sigma-k-lattice");
  FUNC_LIB::function<std::complex<double>, dmn_3<nu, nu, r_HOST>> S_r;  //("Sigma-r-lattice");

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w>> Sigma;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w>> Sigma_stddev;

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w>> Sigma_cluster;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_HOST, w>> Sigma_lattice;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_HOST, w>> Sigma_lattice_interpolated;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_HOST, w>> Sigma_lattice_coarsegrained;

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w>> G_k_w;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w>> G_k_w_stddev;
  FUNC_LIB::function<double, dmn_4<nu, nu, k_DCA, t>> G_k_t;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_DCA, w>> G_r_w;
  FUNC_LIB::function<double, dmn_4<nu, nu, r_DCA, t>> G_r_t;

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w>> G0_k_w;
  FUNC_LIB::function<double, dmn_4<nu, nu, k_DCA, t>> G0_k_t;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_DCA, w>> G0_r_w;
  FUNC_LIB::function<double, dmn_4<nu, nu, r_DCA, t>> G0_r_t;

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w>> G0_k_w_cluster_excluded;
  FUNC_LIB::function<double, dmn_4<nu, nu, k_DCA, t>> G0_k_t_cluster_excluded;
  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_DCA, w>> G0_r_w_cluster_excluded;
  FUNC_LIB::function<double, dmn_4<nu, nu, r_DCA, t>> G0_r_t_cluster_excluded;

  FUNC_LIB::function<std::complex<double>, dmn_8<b, b, b, b, k_DCA, k_DCA, w_VERTEX, w_VERTEX>> G4_k_k_w_w;
  FUNC_LIB::function<std::complex<double>, dmn_8<b, b, b, b, k_DCA, k_DCA, w_VERTEX, w_VERTEX>>
      G4_k_k_w_w_stddev;

  FUNC_LIB::function<double, dmn_4<nu, nu, r_DCA, t>> K_r_t;

  FUNC_LIB::function<double, nu> orbital_occupancy;
};

template <class parameters_type>
DCA_data<parameters_type>::DCA_data(parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      H_symmetry("H_symmetry"),
      H_interactions("interaction-matrix"),

      H_DCA("H_DCA"),
      H_HOST("H_HOST"),
      H_LDA("H_LDA"),

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

      // visited_expansion_order_k("<k>"),

      orbital_occupancy("orbital_occupancy")  //,

// mu("chemical_potential_mu")
{
  H_symmetry = -1;
}

template <class parameters_type>
void DCA_data<parameters_type>::read(std::string filename) {
  if (concurrency.id() == 0)
    std::cout << "\n\n\t starts reading \n\n";

  if (concurrency.id() == concurrency.first()) {
    IO::FORMAT FORMAT = parameters.get_output_format();

    switch (FORMAT) {
      case IO::JSON: {
        IO::reader<IO::JSON> reader;

        reader.open_file(filename);

        this->read(reader);

        reader.close_file();
      } break;

      case IO::HDF5: {
        IO::reader<IO::HDF5> reader;

        reader.open_file(filename);

        this->read(reader);

        reader.close_file();
      } break;

      default:
        throw std::logic_error(__FUNCTION__);
    }
  }

  concurrency.broadcast(parameters.get_chemical_potential());

  concurrency.broadcast_object(Sigma);
  // to do: Interpolate Sigma if \beta from file != \beta ?
  if (parameters.do_CPE())
    concurrency.broadcast_object(G_k_t);

  if (parameters.get_vertex_measurement_type() != NONE) {
    concurrency.broadcast_object(G_k_w);
    concurrency.broadcast_object(G4_k_k_w_w);
  }
}

template <class parameters_type>
template <IO::FORMAT DATA_FORMAT>
void DCA_data<parameters_type>::read(IO::reader<DATA_FORMAT>& reader) {
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

    if (parameters.do_CPE())
      reader.execute(G_k_t);

    if (vertex_measurement != "NONE") {
      reader.execute(G_k_w);

      reader.execute(G4_k_k_w_w);
    }

    reader.close_group();
  }
}

template <class parameters_type>
void DCA_data<parameters_type>::write(std::string file_name) {
  IO::FORMAT FORMAT = parameters.get_output_format();

  std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

  switch (FORMAT) {
    case IO::JSON: {
      IO::writer<IO::JSON> writer;
      {
        writer.open_file(file_name);

        parameters.write(writer);
        this->write(writer);

        writer.close_file();
      }
    } break;

    case IO::HDF5: {
      IO::writer<IO::HDF5> writer;
      {
        writer.open_file(file_name);

        parameters.write(writer);
        this->write(writer);

        writer.close_file();
      }
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <class parameters_type>
template <IO::FORMAT DATA_FORMAT>
void DCA_data<parameters_type>::write(IO::writer<DATA_FORMAT>& writer) {
  writer.open_group("functions");

  writer.execute(band_structure);

  if (parameters.do_DCA_plus()) {
    writer.execute(Sigma_band_structure);
    writer.execute(Sigma_cluster_band_structure);
    writer.execute(Sigma_lattice_band_structure);
    writer.execute(Sigma_band_structure_interpolated);
    writer.execute(Sigma_band_structure_coarsegrained);

    writer.execute(S_k);
    writer.execute(S_r);

    writer.execute(G_k);
  }

  if (!parameters.do_DCA_plus()) {  // Compute Sigma-r-DCA for the lowest frequency
                                    // via Fourier transformation of DCA cluster
                                    // Sigma.
    FUNC_LIB::function<std::complex<double>, dmn_3<nu, nu, r_DCA>> S_r_DCA("Sigma-r-DCA");

    FUNC_LIB::function<std::complex<double>, dmn_3<nu, nu, k_DCA>> S_k_DCA("Sigma-k-DCA");
    std::memcpy(&S_k_DCA(0), &Sigma(0, 0, 0, w::dmn_size() / 2),
                sizeof(std::complex<double>) * std::pow(2 * b::dmn_size(), 2.) * k_DCA::dmn_size());
    math_algorithms::functional_transforms::TRANSFORM<k_DCA, r_DCA>::execute(S_k_DCA, S_r_DCA);

    writer.execute(S_r_DCA);
  }

  writer.execute(Sigma);
  writer.execute(Sigma_stddev);

  if (parameters.dump_lattice_Self_energy()) {
    writer.execute(Sigma_lattice);
  }

  if (parameters.do_CPE() or parameters.do_equal_time_measurements() or
      parameters.dump_cluster_Greens_functions()) {
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

  if (parameters.get_vertex_measurement_type() != NONE) {
    if (not(parameters.do_CPE() or parameters.do_equal_time_measurements() or
            parameters.dump_cluster_Greens_functions())) {
      writer.execute(G_k_w);
      writer.execute(G_k_w_stddev);
    }

    writer.execute(G4_k_k_w_w);
    writer.execute(G4_k_k_w_w_stddev);
  }

  writer.close_group();
}

template <class parameters_type>
void DCA_data<parameters_type>::initialize() {
  initialize_H_0_and_H_i();

  initialize_G0();
}

template <class parameters_type>
void DCA_data<parameters_type>::initialize_H_0_and_H_i() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t initialize H_0(k) and H_i " << dca::util::print_time() << "\n";

  parameters_type::model_type::initialize_H_LDA(H_LDA, parameters);

  parameters_type::model_type::initialize_H_interaction(H_interactions, parameters);

  parameters_type::model_type::initialize_H_symmetries(H_symmetry);

  {
    wannier_interpolation<k_LDA, k_DCA>::execute(H_LDA, H_DCA);
    wannier_interpolation<k_LDA, k_HOST>::execute(H_LDA, H_HOST);

    compute_band_structure::execute(parameters, H_LDA, band_structure);
  }

  if (concurrency.id() == concurrency.first())
    std::cout << "\t finished H_0(k) and H_i " << dca::util::print_time() << "\n";
}

template <class parameters_type>
void DCA_data<parameters_type>::initialize_G0() {
  profiler_type prof("initialize-G0", "input", __LINE__);

  if (concurrency.id() == 0)
    std::cout << "\n\n\t initialize G0 " << dca::util::print_time() << "\n";

  DCA::coarsegraining_sp<parameters_type, k_DCA> coarsegrain_obj(parameters);

  if (concurrency.id() == 0)
    std::cout << "\t\t start coarsegraining G0_k_w " << dca::util::print_time() << "\n";

  {
    FUNC_LIB::function<std::complex<double>, nu_nu_k_DCA_w> Sigma_zero;
    Sigma_zero = 0.;

    coarsegrain_obj.compute_G_K_w(H_HOST, Sigma_zero, G0_k_w);

    symmetrize::execute(G0_k_w, H_symmetry, true);
  }

  if (concurrency.id() == 0)
    std::cout << "\t\t start coarsegraining G0_k_t " << dca::util::print_time() << "\n";

  {
    coarsegrain_obj.compute_G0_K_t(H_HOST, G0_k_t);

    symmetrize::execute(G0_k_t, H_symmetry, true);
  }
  // test_initialize_G0();

  if (concurrency.id() == 0)
    std::cout << "\n\t\t FT G0_k_w, G0_k_t --> G0_r_w, G0_r_t " << dca::util::print_time() << "\n";

  {
    math_algorithms::functional_transforms::TRANSFORM<k_DCA, r_DCA>::execute(G0_k_w, G0_r_w);
    math_algorithms::functional_transforms::TRANSFORM<k_DCA, r_DCA>::execute(G0_k_t, G0_r_t);

    symmetrize::execute(G0_r_t, H_symmetry, true);
  }

  if (concurrency.id() == 0)
    std::cout << "\t finished G0 " << dca::util::print_time();
}

template <class parameters_type>
void DCA_data<parameters_type>::initialize_Sigma() {
  profiler_type prof("initialize-Sigma", "input", __LINE__);

  if (parameters.get_Sigma_file() != "zero")
    this->read(parameters.get_Sigma_file());
}

template <class parameters_type>
void DCA_data<parameters_type>::compute_single_particle_properties() {
  {
    std::memcpy(&S_k(0), &Sigma_lattice(0, 0, 0, w::dmn_size() / 2),
                sizeof(std::complex<double>) * std::pow(2 * b::dmn_size(), 2.) * k_HOST::dmn_size());

    math_algorithms::functional_transforms::TRANSFORM<k_HOST, r_HOST>::execute(S_k, S_r);
  }

  {
    int w_ind = w::dmn_size() / 2;

    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> I_k("I_matrix", nu::dmn_size());
    dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_inv("G_inv", nu::dmn_size());

    LIN_ALG::GEINV<dca::linalg::CPU>::plan<std::complex<double>> geinv_obj(G_inv);

    std::complex<double> i_wm_plus_mu;

    i_wm_plus_mu.real(parameters.get_chemical_potential());
    i_wm_plus_mu.imag(w::get_elements()[w_ind]);

    for (int i = 0; i < nu::dmn_size(); i++)
      I_k(i, i) = i_wm_plus_mu;

    for (int k_ind = 0; k_ind < k_HOST::dmn_size(); k_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_inv(i, j) = I_k(i, j) - H_HOST(i, j, k_ind) - Sigma_lattice(i, j, k_ind, w_ind);

      geinv_obj.execute(G_inv);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_k(i, j, k_ind) = G_inv(i, j);
    }
  }
}

template <class parameters_type>
void DCA_data<parameters_type>::compute_Sigma_bands() {
  {
    Sigma_band_structure.reset();
    Sigma_cluster_band_structure.reset();

    std::vector<std::pair<double, int>> length_and_distance(k_DCA::dmn_size(),
                                                            std::pair<double, int>(0, -1));

    for (int k_ind = 0; k_ind < k_domain_cut_dmn_type::dmn_size(); ++k_ind) {
      std::vector<double> k_vec = cluster_operations::translate_inside_cluster(
          k_domain_cut_dmn_type::get_elements()[k_ind],
          DCA_k_cluster_type::get_super_basis_vectors());

      for (int K_ind = 0; K_ind < k_DCA::dmn_size(); ++K_ind) {
        length_and_distance[K_ind].second = K_ind;

        length_and_distance[K_ind].first = cluster_operations::minimal_distance(
            k_vec, k_DCA::get_elements()[K_ind], DCA_k_cluster_type::get_super_basis_vectors());
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
  if (parameters.do_DCA_plus()) {
    FUNC_LIB::function<std::complex<double>, dmn_2<nu, k_HOST>> S_k_dmn("S_k_dmn_s");

    for (int b_ind = 0; b_ind < b::dmn_size(); ++b_ind)
      for (int s_ind = 0; s_ind < s::dmn_size(); ++s_ind)
        for (int k_ind = 0; k_ind < k_HOST::dmn_size(); ++k_ind)
          S_k_dmn(b_ind, s_ind, k_ind) =
              Sigma_lattice(b_ind, s_ind, b_ind, s_ind, k_ind, w::dmn_size() / 2);

    hspline_interpolation<k_HOST, k_domain_cut_dmn_type>::execute(
        S_k_dmn, Sigma_lattice_band_structure, -1. / 2.);
  }

  Sigma_band_structure_interpolated.reset();
  if (true) {
    FUNC_LIB::function<std::complex<double>, dmn_2<nu, k_HOST>> S_k_dmn("S_k_dmn_s");

    for (int b_ind = 0; b_ind < b::dmn_size(); ++b_ind)
      for (int s_ind = 0; s_ind < s::dmn_size(); ++s_ind)
        for (int k_ind = 0; k_ind < k_HOST::dmn_size(); ++k_ind)
          S_k_dmn(b_ind, s_ind, k_ind) =
              Sigma_lattice_interpolated(b_ind, s_ind, b_ind, s_ind, k_ind, w::dmn_size() / 2);

    hspline_interpolation<k_HOST, k_domain_cut_dmn_type>::execute(
        S_k_dmn, Sigma_band_structure_interpolated, -1. / 2.);
  }

  Sigma_band_structure_coarsegrained.reset();
  if (parameters.do_DCA_plus()) {
    FUNC_LIB::function<std::complex<double>, dmn_2<nu, k_HOST>> S_k_dmn("S_k_dmn_s");

    for (int b_ind = 0; b_ind < b::dmn_size(); ++b_ind)
      for (int s_ind = 0; s_ind < s::dmn_size(); ++s_ind)
        for (int k_ind = 0; k_ind < k_HOST::dmn_size(); ++k_ind)
          S_k_dmn(b_ind, s_ind, k_ind) =
              Sigma_lattice_coarsegrained(b_ind, s_ind, b_ind, s_ind, k_ind, w::dmn_size() / 2);

    hspline_interpolation<k_HOST, k_domain_cut_dmn_type>::execute(
        S_k_dmn, Sigma_band_structure_coarsegrained, -1. / 2.);
  }
}

template <class parameters_type>
void DCA_data<parameters_type>::print_Sigma_QMC_versus_Sigma_cg() {
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

    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); ++k_ind) {
      VECTOR_OPERATIONS::PRINT(k_DCA::get_elements()[k_ind]);
      std::cout << real(Sigma(0, 0, k_ind, w::dmn_size() / 2)) << "\t"
                << imag(Sigma(0, 0, k_ind, w::dmn_size() / 2)) << "\t";
      std::cout << real(Sigma_cluster(0, 0, k_ind, w::dmn_size() / 2)) << "\t"
                << imag(Sigma_cluster(0, 0, k_ind, w::dmn_size() / 2)) << "\n";
    }
    std::cout << "\n\n";
  }
}
}

#endif  // PHYS_LIBRARY_DCA_DATA_DCA_DATA_H
