// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes susceptibilities by using the BseClusterSolver and BseLatticeSolver classes.

#ifndef DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_SOLVER_HPP
#define DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_SOLVER_HPP

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_writer.hpp"
#include "dca/phys/dca_algorithms/compute_band_structure.hpp"
#include "dca/phys/dca_analysis/bse_solver/bse_cluster_solver.hpp"
#include "dca/phys/dca_analysis/bse_solver/bse_lattice_solver.hpp"
#include "dca/phys/dca_step/symmetrization/diagrammatic_symmetries.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/brillouin_zone_path_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

template <typename ParametersType, typename DcaDataType>
class BseSolver {
public:
  using scalartype = double;

  using profiler_t = typename ParametersType::profiler_type;
  using concurrency_t = typename ParametersType::concurrency_type;

  const static int N_LAMBDAS = 10;
  using lambda_dmn_type = func::dmn_0<func::dmn<N_LAMBDAS, int>>;

  const static int N_HARMONICS = 3;
  using harmonics_dmn_type = func::dmn_0<func::dmn<N_HARMONICS, int>>;

  using w_VERTEX = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using b_b = func::dmn_variadic<b, b>;

  using k_DCA =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::CLUSTER,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  using k_LDA =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::LATTICE_SP,
                                          domains::MOMENTUM_SPACE, domains::PARALLELLEPIPEDUM>>;
  using k_HOST =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::LATTICE_SP,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  using k_HOST_VERTEX =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::LATTICE_TP,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  using k_dmn_cut_type = func::dmn_0<domains::brillouin_zone_path_domain<domains::SQUARE_2D_LATTICE>>;

  using cluster_eigenvector_dmn_t = func::dmn_variadic<b, b, k_DCA, w_VERTEX>;
  using lattice_eigenvector_dmn_t = func::dmn_variadic<b, b, k_HOST_VERTEX, w_VERTEX>;

  using DCA_matrix_dmn_t = func::dmn_variadic<func::dmn_variadic<b, b, k_DCA, w_VERTEX>,
                                              func::dmn_variadic<b, b, k_DCA, w_VERTEX>>;
  using HOST_matrix_dmn_t = func::dmn_variadic<func::dmn_variadic<b, b, k_HOST_VERTEX, w_VERTEX>,
                                               func::dmn_variadic<b, b, k_HOST_VERTEX, w_VERTEX>>;

  BseSolver(ParametersType& parameters, DcaDataType& MOMS);

  void write();
  template <typename Writer>
  void write(Writer& write);

  void calculate_susceptibilities_2();

  func::function<std::complex<scalartype>, lambda_dmn_type>& get_leading_eigenvalues() {
    return BSE_lattice_solver_obj.get_leading_eigenvalues();
  };

private:
  void initialize_wave_functions();

  void apply_symmetries();

  ParametersType& parameters;
  concurrency_t& concurrency;

  DcaDataType& MOMS;

  BseClusterSolver<ParametersType, DcaDataType> BSE_cluster_solver_obj;
  BseLatticeSolver<ParametersType, DcaDataType> BSE_lattice_solver_obj;

  cluster_eigenvector_dmn_t cluster_eigenvector_dmn;
  lattice_eigenvector_dmn_t lattice_eigenvector_dmn;

  diagrammatic_symmetries<ParametersType> diagrammatic_symmetries_obj;

  func::function<std::string, harmonics_dmn_type> wave_functions_names;
  func::function<std::complex<scalartype>, func::dmn_variadic<k_HOST_VERTEX, harmonics_dmn_type>> harmonics;

  func::function<std::complex<double>, DCA_matrix_dmn_t> G4;
  func::function<std::complex<double>, DCA_matrix_dmn_t> G4_0;

  func::function<std::complex<scalartype>, DCA_matrix_dmn_t> Gamma_cluster;
  func::function<std::complex<scalartype>, HOST_matrix_dmn_t> Gamma_lattice;

  func::function<std::complex<double>, func::dmn_variadic<b_b, b_b, k_HOST_VERTEX, w_VERTEX>>
      chi_0_function;  //("phi");

  func::function<std::complex<scalartype>, HOST_matrix_dmn_t> chi;
  func::function<std::complex<scalartype>, HOST_matrix_dmn_t> chi_0;

  func::function<std::complex<scalartype>, func::dmn_0<func::dmn<1, int>>> chi_q;

  func::function<std::complex<scalartype>, harmonics_dmn_type> P_q_cluster;
  func::function<std::complex<scalartype>, harmonics_dmn_type> P_q_lattice;

  func::function<std::complex<scalartype>, lambda_dmn_type> leading_eigenvalues;
  func::function<std::complex<scalartype>, lambda_dmn_type> leading_phi_t_chi_0_phi;
  func::function<std::complex<scalartype>, lambda_dmn_type> leading_phi_t_Gamma_phi;

  func::function<std::complex<scalartype>, func::dmn_variadic<lambda_dmn_type, harmonics_dmn_type>>
      leading_symmetries;
  func::function<std::complex<scalartype>, func::dmn_variadic<lambda_dmn_type, lattice_eigenvector_dmn_t>>
      leading_eigenvectors;

  func::function<std::complex<scalartype>, k_dmn_cut_type> S_k_cut;
  func::function<std::complex<scalartype>, k_dmn_cut_type> a_k_cut;

  func::function<std::complex<scalartype>, func::dmn_variadic<lambda_dmn_type, cluster_eigenvector_dmn_t>>
      leading_U_K;
  func::function<std::complex<scalartype>, func::dmn_variadic<lambda_dmn_type, cluster_eigenvector_dmn_t>>
      leading_Vt_K;

  func::function<std::complex<scalartype>, func::dmn_variadic<lambda_dmn_type, lattice_eigenvector_dmn_t>>
      leading_U_k;
  func::function<std::complex<scalartype>, func::dmn_variadic<lambda_dmn_type, lattice_eigenvector_dmn_t>>
      leading_Vt_k;
};

template <typename ParametersType, typename DcaDataType>
BseSolver<ParametersType, DcaDataType>::BseSolver(ParametersType& parameters_ref,
                                                  DcaDataType& MOMS_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      MOMS(MOMS_ref),

      BSE_cluster_solver_obj(parameters, MOMS),
      BSE_lattice_solver_obj(parameters, MOMS),

      cluster_eigenvector_dmn(),
      lattice_eigenvector_dmn(),

      diagrammatic_symmetries_obj(parameters),

      harmonics("harmonics"),

      G4("G4"),
      G4_0("G4_0"),

      Gamma_cluster("Gamma_cluster"),
      Gamma_lattice("Gamma_lattice"),

      chi_0_function("chi_0_function"),

      chi("chi"),
      chi_0("chi_0"),

      chi_q("chi_q"),

      P_q_cluster("P_q_cluster"),
      P_q_lattice("P_q_lattice"),

      leading_eigenvalues("leading_eigenvalues"),
      leading_phi_t_chi_0_phi("leading_phi_t_chi_0_phi"),
      leading_phi_t_Gamma_phi("leading_phi_t_Gamma_phi"),

      leading_symmetries("leading_symmetries"),
      leading_eigenvectors("leading_eigenvectors"),

      S_k_cut("S_k_cut"),
      a_k_cut("a_k_cut"),

      leading_U_K("leading_U_K"),
      leading_Vt_K("leading_Vt_K"),

      leading_U_k("leading_U_k_interpolated"),
      leading_Vt_k("leading_Vt_k_interpolated") {
  initialize_wave_functions();

  {
    profiler_t prof("compute-H(k)", "input", __LINE__);

    ParametersType::model_type::initialize_H_0(parameters, MOMS.H_DCA);
    ParametersType::model_type::initialize_H_0(parameters, MOMS.H_HOST);
  }

  {
    profiler_t prof("compute-band-structure", "input", __LINE__);
    compute_band_structure::execute(parameters, MOMS.band_structure);
  }
}

template <typename ParametersType, typename DcaDataType>
void BseSolver<ParametersType, DcaDataType>::write() {
  const std::string& output_format = parameters.get_output_format();
  const std::string& file_name =
      parameters.get_directory() + parameters.get_susceptibilities_file_name();

  std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

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
    // MOMS.write(writer);
    this->write(writer);

    writer.close_file();
  }

  else
    throw std::logic_error(__FUNCTION__);
}

template <typename ParametersType, typename DcaDataType>
template <typename Writer>
void BseSolver<ParametersType, DcaDataType>::write(Writer& writer) {
  writer.open_group("analysis-functions");

  {
    BSE_lattice_solver_obj.write(writer);
    BSE_cluster_solver_obj.write(writer);
  }

  writer.close_group();
}

template <typename ParametersType, typename DcaDataType>
void BseSolver<ParametersType, DcaDataType>::initialize_wave_functions() {
  wave_functions_names = "no-name";

  wave_functions_names(0) = "s-wave";
  wave_functions_names(1) = "p-wave";
  wave_functions_names(2) = "d-wave";

  {  // s-wave
    std::complex<scalartype> norm_psi = 0;

    for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
      harmonics(k_ind, 0) = 1.;

    for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
      norm_psi += harmonics(k_ind, 0) * conj(harmonics(k_ind, 0));

    for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
      harmonics(k_ind, 0) /= std::sqrt(norm_psi);
  }

  {  // p-wave
    std::complex<scalartype> norm_psi = 0;

    scalartype alpha_x = 1;  // host_vertex_cluster_type::get_r_basis()[0][0];
    scalartype alpha_y = 1;  // host_vertex_cluster_type::get_r_basis()[1][1];

    for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
      harmonics(k_ind, 1) = (cos(alpha_x * k_HOST_VERTEX::get_elements()[k_ind][0]) +
                             cos(alpha_y * k_HOST_VERTEX::get_elements()[k_ind][1]));

    for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
      norm_psi += harmonics(k_ind, 1) * conj(harmonics(k_ind, 1));

    for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
      harmonics(k_ind, 1) /= std::sqrt(norm_psi);
  }

  {  // d-wave
    std::complex<scalartype> norm_psi = 0;

    scalartype alpha_x = 1;  // host_vertex_cluster_type::get_r_basis()[0][0];
    scalartype alpha_y = 1;  // host_vertex_cluster_type::get_r_basis()[1][1];

    for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
      harmonics(k_ind, 2) = (cos(alpha_x * k_HOST_VERTEX::get_elements()[k_ind][0]) -
                             cos(alpha_y * k_HOST_VERTEX::get_elements()[k_ind][1]));

    for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
      norm_psi += harmonics(k_ind, 2) * conj(harmonics(k_ind, 2));

    for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
      harmonics(k_ind, 2) /= std::sqrt(norm_psi);
  }
}

template <typename ParametersType, typename DcaDataType>
void BseSolver<ParametersType, DcaDataType>::calculate_susceptibilities_2() {
  if (concurrency.id() == concurrency.last())
    std::cout << "\t" << __FUNCTION__ << std::endl;

  BSE_cluster_solver_obj.compute_Gamma_cluster();

  Gamma_cluster = BSE_cluster_solver_obj.get_Gamma_matrix();

  BSE_lattice_solver_obj.compute_Gamma_lattice_3(Gamma_cluster);
  BSE_lattice_solver_obj.compute_chi_0_lattice(chi_0);

  Gamma_lattice = BSE_lattice_solver_obj.get_Gamma_lattice();

  BSE_lattice_solver_obj.diagonalize_Gamma_chi_0(Gamma_lattice, chi_0);
}

template <typename ParametersType, typename DcaDataType>
void BseSolver<ParametersType, DcaDataType>::apply_symmetries() {
  profiler_t prof(__FUNCTION__, __FILE__, __LINE__);

  if (concurrency.id() == concurrency.last())
    std::cout << "\t" << __FUNCTION__ << std::endl << std::endl;

  symmetrize::execute(MOMS.Sigma, MOMS.H_symmetry);

  symmetrize::execute(MOMS.G_k_w, MOMS.H_symmetry);
}

}  // analysis
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_SOLVER_HPP
