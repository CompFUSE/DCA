// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides an exact diagonalization (ED) cluster solver.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ED_CLUSTER_SOLVER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ED_CLUSTER_SOLVER_HPP

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_writer.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_data/dca_data_real_freq.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/fermionic_overlap_matrices.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/fock_space.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/greens_functions/sp_greens_function.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/greens_functions/tp_greens_function.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/hamiltonian.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/options.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_real_axis.hpp"
#include "dca/phys/four_point_type.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
class EDClusterSolver {
public:
  using MOMS_w_imag_type = DcaData<parameters_type>;
  using MOMS_w_real_type = DcaDataRealFreq<parameters_type>;

  using ed_options_type = ed::Options<parameters_type>;

  using b = typename ed_options_type::b;
  using nu = typename ed_options_type::nu;

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  using w_REAL = func::dmn_0<domains::frequency_domain_real_axis>;

public:
  EDClusterSolver(parameters_type& parameters_ref, MOMS_type& MOMS_ref,
                  MOMS_w_real_type& MOMS_real_ref);

  void initialize(int dca_iteration);

  void execute();

  template <typename dca_info_struct_t>
  void finalize(dca_info_struct_t& dca_info_struct);

  void write(std::string file_name);

private:
  parameters_type& parameters;
  typename parameters_type::concurrency_type& concurrency;
  MOMS_type& MOMS_imag;
  MOMS_w_real_type& MOMS_real;

  ed::Fock_space<parameters_type, ed_options_type> Fock_obj;
  ed::Hamiltonian<parameters_type, ed_options_type> Ham_obj;
  ed::fermionic_overlap_matrices<parameters_type, ed_options_type> overlap_obj;
  ed::SpGreensFunction<parameters_type, ed_options_type> sp_Greens_function_obj;
  ed::TpGreensFunction<parameters_type, ed_options_type> tp_Greens_function_obj;
};

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
EDClusterSolver<device_t, parameters_type, MOMS_type>::EDClusterSolver(
    parameters_type& parameters_ref, MOMS_type& MOMS_imag_ref, MOMS_w_real_type& MOMS_real_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      MOMS_imag(MOMS_imag_ref),
      MOMS_real(MOMS_real_ref),

      Fock_obj(true, true),
      Ham_obj(parameters),
      overlap_obj(parameters, Ham_obj),

      sp_Greens_function_obj(parameters, Ham_obj, overlap_obj),
      tp_Greens_function_obj(parameters, Ham_obj, overlap_obj) {
#ifdef DEBUG
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n\nFock-space without symmetries:" << std::endl;
    Fock_obj.print_subspaces();
  }
#endif

  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n\n\n"
              << "Apply translation symmetry ..." << std::endl;
  }
  Fock_obj.apply_translation_symmetry();

  if (concurrency.id() == concurrency.first()) {
    std::cout << dca::util::print_time() << std::endl;
    std::cout << "\nCreate representation ..." << std::endl;
  }
  Fock_obj.initialize_rep();

  if (concurrency.id() == concurrency.first()) {
    std::cout << dca::util::print_time() << std::endl;
  }

#ifdef DEBUG
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\nFock-space with symmetries:" << std::endl;
    Fock_obj.print_subspaces();
  }
#endif
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void EDClusterSolver<device_t, parameters_type, MOMS_type>::initialize(int /*dca_iteration*/) {
  func::function<std::complex<double>, func::dmn_variadic<nu, nu, RClusterDmn>> H_DCA;

  math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(MOMS_imag.H_DCA, H_DCA);

  Ham_obj.initialize(H_DCA, MOMS_imag.H_interactions);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void EDClusterSolver<device_t, parameters_type, MOMS_type>::execute() {
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n" << __FUNCTION__ << "\n" << std::endl;
  }

  // creation and annihilation matrices
  overlap_obj.construct_creation_set_all();
  overlap_obj.construct_annihilation_set_all();
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n" << dca::util::print_time() << std::endl;
  }
  overlap_obj.construct_creation_set_nonzero_sparse();
  overlap_obj.construct_annihilation_set_nonzero_sparse();

  {  // non-interacting Greensfunction
    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n" << dca::util::print_time() << std::endl;
    }
    Ham_obj.construct_Hamiltonians(false);

    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n" << dca::util::print_time() << std::endl;
    }
    Ham_obj.diagonalize_Hamiltonians_st();

    if (concurrency.id() == concurrency.first()) {
      std::cout << dca::util::print_time() << std::endl;
    }
    Ham_obj.set_spectrum(MOMS_real.E0_w);

    sp_Greens_function_obj.compute_all_sp_functions_slow(MOMS_imag, MOMS_real, false);

    if (parameters.get_four_point_type() != NONE) {
      if (concurrency.id() == concurrency.first()) {
        std::cout << "\n" << dca::util::print_time() << "\n" << std::endl;
      }
      tp_Greens_function_obj.compute_two_particle_Greens_function(false);
      if (concurrency.id() == concurrency.first()) {
        std::cout << "\n" << dca::util::print_time() << "\n" << std::endl;
      }
      tp_Greens_function_obj.compute_particle_particle_superconducting_A(MOMS_imag.get_G4());
    }
  }

  // interacting Greensfunction
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n" << dca::util::print_time() << std::endl;
  }
  Ham_obj.construct_Hamiltonians(true);

  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n" << dca::util::print_time() << std::endl;
  }
  Ham_obj.diagonalize_Hamiltonians_st();

  if (concurrency.id() == concurrency.first()) {
    std::cout << dca::util::print_time() << std::endl;
  }
  Ham_obj.set_spectrum(MOMS_real.E_w);

  sp_Greens_function_obj.compute_all_sp_functions_slow(MOMS_imag, MOMS_real, true);

  if (parameters.get_four_point_type() != NONE) {
    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n" << dca::util::print_time() << "\n" << std::endl;
    }
    tp_Greens_function_obj.compute_two_particle_Greens_function(true);
    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n" << dca::util::print_time() << "\n" << std::endl;
    }
    tp_Greens_function_obj.compute_particle_particle_superconducting_A(MOMS_imag.get_G4());
  }

  MOMS_real.A_w = 0;
  for (int l = 0; l < w_REAL::dmn_size(); l++)
    for (int j = 0; j < KClusterDmn::dmn_size(); j++)
      for (int i = 0; i < 2 * b::dmn_size(); i++)
        MOMS_real.A_w(l) -= imag(MOMS_real.G_k_w(i, i, j, l));

  MOMS_real.A0_w = 0;
  for (int l = 0; l < w_REAL::dmn_size(); l++)
    for (int j = 0; j < KClusterDmn::dmn_size(); j++)
      for (int i = 0; i < 2 * b::dmn_size(); i++)
        MOMS_real.A0_w(l) -= imag(MOMS_real.G0_k_w(i, i, j, l));

  MOMS_real.A_w *= 1. / double(M_PI * KClusterDmn::dmn_size() * 2 * b::dmn_size());
  MOMS_real.A0_w *= 1. / double(M_PI * KClusterDmn::dmn_size() * 2 * b::dmn_size());

  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n" << dca::util::print_time() << std::endl;
  }

  sp_Greens_function_obj.compute_S_k_w(MOMS_imag.G_k_w, MOMS_imag.G0_k_w, MOMS_imag.Sigma);
  sp_Greens_function_obj.compute_S_k_w(MOMS_real.G_k_w, MOMS_real.G0_k_w, MOMS_real.Sigma);

  if (concurrency.id() == concurrency.first()) {
    std::cout << dca::util::print_time() << std::endl;
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
template <typename dca_info_struct_t>
void EDClusterSolver<device_t, parameters_type, MOMS_type>::finalize(
    dca_info_struct_t& /*dca_info_struct*/) {
  for (int l = 0; l < MOMS_imag.G_r_w.size(); l++)
    MOMS_imag.Sigma_cluster(l) = MOMS_imag.Sigma(l);

  for (int l = 0; l < MOMS_imag.G0_k_w.size(); l++)
    MOMS_imag.G0_k_w_cluster_excluded(l) = MOMS_imag.G0_k_w(l);

  for (int l = 0; l < MOMS_imag.G0_r_w.size(); l++)
    MOMS_imag.G0_r_w_cluster_excluded(l) = MOMS_imag.G0_r_w(l);

  for (int l = 0; l < MOMS_imag.G0_k_t.size(); l++)
    MOMS_imag.G0_k_t_cluster_excluded(l) = MOMS_imag.G0_k_t(l);

  for (int l = 0; l < MOMS_imag.G0_r_t.size(); l++)
    MOMS_imag.G0_r_t_cluster_excluded(l) = MOMS_imag.G0_r_t(l);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void EDClusterSolver<device_t, parameters_type, MOMS_type>::write(std::string file_name) {
  std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

  const std::string& output_format = parameters.get_output_format();

  if (output_format == "JSON") {
    dca::io::JSONWriter writer;
    writer.open_file(file_name);

    parameters.write(writer);
    MOMS_imag.write(writer);
    MOMS_real.write(writer);

    if (parameters.get_four_point_type() != NONE) {
      std::cout << "\n\n\t\t start writing tp-Greens-function\n\n";
      tp_Greens_function_obj.write(writer);
    }

    writer.close_file();
  }

  else if (output_format == "HDF5") {
    dca::io::HDF5Writer writer;
    writer.open_file(file_name);

    std::cout << "\n\n\t\t start writing parameters\n\n";
    parameters.write(writer);

    std::cout << "\n\n\t\t start writing MOMS_imag\n\n";
    MOMS_imag.write(writer);

    std::cout << "\n\n\t\t start writing MOMS_real\n\n";
    MOMS_real.write(writer);

    if (parameters.get_four_point_type() != NONE) {
      std::cout << "\n\n\t\t start writing tp-Greens-function\n\n";
      tp_Greens_function_obj.write(writer);
    }

    writer.close_file();
  }

  else
    throw std::logic_error(__FUNCTION__);
}

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ED_CLUSTER_SOLVER_HPP
