// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// This class computes susceptibilities for inelastic neutron scattering.
//
// TODO: Add descriptions to (public) methods.

#ifndef DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_SOLVER_EXT_HPP
#define DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_SOLVER_EXT_HPP

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
#include "dca/phys/dca_analysis/bse_solver/bse_cluster_solver_ext.hpp"
#include "dca/phys/dca_analysis/bse_solver/bse_lattice_solver_ext.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

template <typename ParametersType, typename DcaDataType>
class BseSolverExt {
public:
  using CDA = ClusterDomainAliases<ParametersType::lattice_type::DIMENSION>;
  using Scalar = typename ParametersType::Scalar;
  using Real = typename ParametersType::Real;
  using Complex = typename std::complex<typename ParametersType::Real>;

  using ProfilerType = typename ParametersType::profiler_type;
  using ConcurrencyType = typename ParametersType::concurrency_type;

  using BseClusterSolverType = BseClusterSolverExt<ParametersType, DcaDataType, Scalar>;
  using BseLatticeSolverType = BseLatticeSolverExt<ParametersType, DcaDataType, Scalar>;
  using LeadingEigDmn = typename BseLatticeSolverType::LeadingEigDmn;
  using LatticeEigenvectorDmn = typename BseLatticeSolverType::LatticeEigenvectorDmn;

  static constexpr int num_harmonics = 3;
  using HarmonicsDmn = func::dmn_0<func::dmn<num_harmonics, int>>;

  using TpHostKDmn =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::LATTICE_TP,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;

  BseSolverExt(ParametersType& parameters, DcaDataType& dca_data);

  void write();
  template <typename Writer>
  void write(Writer& write);

  void calculateSusceptibilities();

  auto& get_leading_eigenvalues() /*const*/ {
    return bse_lattice_solver_.get_leading_eigenvalues();
  };
  auto& get_leading_eigenvectors() /*const*/ {
    return bse_lattice_solver_.get_leading_eigenvectors();
  };
  auto& get_leading_symmetry_decomposition() /*const*/ {
    return bse_lattice_solver_.get_leading_symmetry_decomposition();
  };

private:
  // TODO: Cleanup.
  void initialize_wave_functions();

  ParametersType& parameters_;
  ConcurrencyType& concurrency_;

  DcaDataType& dca_data_;

  BseClusterSolverType bse_cluster_solver_;
  BseLatticeSolverType bse_lattice_solver_;

  func::function<std::string, HarmonicsDmn> wave_functions_names_;
  func::function<std::complex<Scalar>, func::dmn_variadic<TpHostKDmn, HarmonicsDmn>> harmonics_;
};

template <typename ParametersType, typename DcaDataType>
BseSolverExt<ParametersType, DcaDataType>::BseSolverExt(ParametersType& parameters,
                                                        DcaDataType& dca_data)
    : parameters_(parameters),
      concurrency_(parameters.get_concurrency()),

      dca_data_(dca_data),

      bse_cluster_solver_(parameters, dca_data),
      bse_lattice_solver_(parameters, dca_data),

      wave_functions_names_("wave-functions-names"),
      harmonics_("harmonics") {
  initialize_wave_functions();

  {
    ProfilerType prof("compute-H_0(k)", "input", __LINE__);
    ParametersType::model_type::initializeH0(parameters_, dca_data_.H_DCA);
    ParametersType::model_type::initializeH0(parameters_, dca_data_.H_HOST);
  }

  {
    ProfilerType prof("compute-band-structure", "input", __LINE__);
    compute_band_structure<ParametersType>::execute(parameters_, dca_data_.band_structure);
  }
}

template <typename ParametersType, typename DcaDataType>
void BseSolverExt<ParametersType, DcaDataType>::write() {
  const std::string& output_format = parameters_.get_output_format();
  const std::string& file_name = parameters_.get_directory() + parameters_.get_filename_analysis();

  std::cout << "Start writing " << file_name << "." << std::endl;

  if (output_format == "JSON") {
    dca::io::JSONWriter writer;
    writer.open_file(file_name);

    parameters_.write(writer);
    // dca_data_.write(writer);
    this->write(writer);

    writer.close_file();
  }
  else if (output_format == "HDF5") {
    dca::io::HDF5Writer writer;
    writer.open_file(file_name);

    parameters_.write(writer);
    // dca_data_.write(writer);
    this->write(writer);

    writer.close_file();
  }
  else
    throw std::logic_error(__FUNCTION__);
}

template <typename ParametersType, typename DcaDataType>
template <typename Writer>
void BseSolverExt<ParametersType, DcaDataType>::write(Writer& writer) {
  writer.open_group("analysis-functions");

  bse_lattice_solver_.write(writer);
  bse_cluster_solver_.write(writer);

  writer.close_group();
}

template <typename ParametersType, typename DcaDataType>
void BseSolverExt<ParametersType, DcaDataType>::initialize_wave_functions() {
  wave_functions_names_ = "no-name";

  wave_functions_names_(0) = "s-wave";
  wave_functions_names_(1) = "p-wave";  // cos(kx)+cos(ky)
  wave_functions_names_(2) = "d-wave";  // cos(kx)-cos(ky)

  {  // s-wave
    std::complex<Scalar> norm_psi = 0;

    for (int k_ind = 0; k_ind < TpHostKDmn::dmn_size(); k_ind++)
      harmonics_(k_ind, 0) = 1.;

    for (int k_ind = 0; k_ind < TpHostKDmn::dmn_size(); k_ind++)
      norm_psi += harmonics_(k_ind, 0) * conj(harmonics_(k_ind, 0));

    for (int k_ind = 0; k_ind < TpHostKDmn::dmn_size(); k_ind++)
      harmonics_(k_ind, 0) /= std::sqrt(norm_psi);
  }

  {  // p-wave
    Complex norm_psi = 0;

    Scalar alpha_x = 1;  // host_vertex_cluster_type::get_r_basis()[0][0];
    Scalar alpha_y = 1;  // host_vertex_cluster_type::get_r_basis()[1][1];

    for (int k_ind = 0; k_ind < TpHostKDmn::dmn_size(); k_ind++)
      harmonics_(k_ind, 1) = (cos(alpha_x * TpHostKDmn::get_elements()[k_ind][0]) +
                              cos(alpha_y * TpHostKDmn::get_elements()[k_ind][1]));

    for (int k_ind = 0; k_ind < TpHostKDmn::dmn_size(); k_ind++)
      norm_psi += harmonics_(k_ind, 1) * conj(harmonics_(k_ind, 1));

    for (int k_ind = 0; k_ind < TpHostKDmn::dmn_size(); k_ind++)
      harmonics_(k_ind, 1) /= std::sqrt(norm_psi);
  }

  {  // d-wave
    Complex norm_psi = 0;

    Scalar alpha_x = 1;  // host_vertex_cluster_type::get_r_basis()[0][0];
    Scalar alpha_y = 1;  // host_vertex_cluster_type::get_r_basis()[1][1];

    for (int k_ind = 0; k_ind < TpHostKDmn::dmn_size(); k_ind++)
      harmonics_(k_ind, 2) = (cos(alpha_x * TpHostKDmn::get_elements()[k_ind][0]) -
                              cos(alpha_y * TpHostKDmn::get_elements()[k_ind][1]));

    for (int k_ind = 0; k_ind < TpHostKDmn::dmn_size(); k_ind++)
      norm_psi += harmonics_(k_ind, 2) * conj(harmonics_(k_ind, 2));

    for (int k_ind = 0; k_ind < TpHostKDmn::dmn_size(); k_ind++)
      harmonics_(k_ind, 2) /= std::sqrt(norm_psi);
  }
}

template <typename ParametersType, typename DcaDataType>
void BseSolverExt<ParametersType, DcaDataType>::calculateSusceptibilities() {
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n" << __FUNCTION__ << std::endl;

  bse_cluster_solver_.compute_Gamma_cluster();

  bse_lattice_solver_.computeGammaLattice(bse_cluster_solver_.get_Gamma_cluster());

  if (CDA::KClusterDmn::dmn_size() == 1) {
    bse_lattice_solver_.computeChi0LatticeOverHost();
    bse_lattice_solver_.computeG4LatticeOverHost();
    bse_lattice_solver_.computeChiDblPrime_q_w();
  }
}

}  // namespace analysis
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_SOLVER_EXT_HPP
