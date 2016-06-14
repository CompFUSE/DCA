// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Description

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_CLUSTER_SOLVER_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_CLUSTER_SOLVER_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_template.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "dca/util/print_time.hpp"
#include "comp_library/function_library/include_function_library.h"
#include "comp_library/IO_library/IO.hpp"
#include "comp_library/linalg/linalg_device_types.h"
#include "phys_library/DCA+_data/DCA_data.h"
#include "phys_library/DCA+_data/moms_w_real.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_Fock_space.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_Hamiltonian.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_Greens_functions/sp_Greens_function_advanced_ed.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_Greens_functions/tp_Greens_function_advanced_ed.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_fermionic_ed_type_definitions.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/overlap_matrix.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_real_axis.h"

namespace DCA {

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
class cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type> {
public:
  using MOMS_w_imag_type = DCA_data<parameters_type>;
  using MOMS_w_real_type = MOMS_w_real<parameters_type>;

  using ed_options_type = ADVANCED_EXACT_DIAGONALIZATION::advanced_ed_options<parameters_type>;

  using b = typename ed_options_type::b;
  using nu = typename ed_options_type::nu;
  using r_DCA = typename ed_options_type::r_DCA;
  using k_DCA = typename ed_options_type::k_DCA;

  using w_REAL = dmn_0<frequency_domain_real_axis>;

public:
  cluster_solver(parameters_type& parameters_ref, MOMS_type& MOMS_ref,
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

  ADVANCED_EXACT_DIAGONALIZATION::Fock_space<parameters_type, ed_options_type> Fock_obj;

  ADVANCED_EXACT_DIAGONALIZATION::fermionic_Hamiltonian<parameters_type, ed_options_type> Ham_obj;

  ADVANCED_EXACT_DIAGONALIZATION::fermionic_overlap_matrices<parameters_type, ed_options_type> overlap_obj;

  ADVANCED_EXACT_DIAGONALIZATION::fermionic_sp_Greens_function<parameters_type, ed_options_type>
      sp_Greens_function_obj;
  ADVANCED_EXACT_DIAGONALIZATION::fermionic_tp_Greens_function<parameters_type, ed_options_type>
      tp_Greens_function_obj;
};

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::cluster_solver(
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
  Fock_obj.apply_translation_symmetry(parameters.get_ED_method());

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

  if (parameters.do_orthogonality_check())
    std::cout << "subspaces orthogonal: " << Fock_obj.check_orthogonality() << std::endl;
}

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
void cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::initialize(
    int /*dca_iteration*/) {
  FUNC_LIB::function<std::complex<double>, dmn_3<nu, nu, r_DCA>> H_DCA;

  math_algorithms::functional_transforms::TRANSFORM<k_DCA, r_DCA>::execute(MOMS_imag.H_DCA, H_DCA);

  Ham_obj.initialize(H_DCA, MOMS_imag.H_interactions);
}

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
void cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::execute() {
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

    if (parameters.get_vertex_measurement_type() != NONE) {
      if (concurrency.id() == concurrency.first()) {
        std::cout << "\n" << dca::util::print_time() << "\n" << std::endl;
      }
      tp_Greens_function_obj.compute_two_particle_Greens_function(false);
      if (concurrency.id() == concurrency.first()) {
        std::cout << "\n" << dca::util::print_time() << "\n" << std::endl;
      }
      tp_Greens_function_obj.compute_particle_particle_superconducting_A(MOMS_imag.G4_k_k_w_w);
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

  if (parameters.get_vertex_measurement_type() != NONE) {
    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n" << dca::util::print_time() << "\n" << std::endl;
    }
    tp_Greens_function_obj.compute_two_particle_Greens_function(true);
    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n" << dca::util::print_time() << "\n" << std::endl;
    }
    tp_Greens_function_obj.compute_particle_particle_superconducting_A(MOMS_imag.G4_k_k_w_w);
  }

  MOMS_real.A_w = 0;
  for (int l = 0; l < w_REAL::dmn_size(); l++)
    for (int j = 0; j < k_DCA::dmn_size(); j++)
      for (int i = 0; i < 2 * b::dmn_size(); i++)
        MOMS_real.A_w(l) -= imag(MOMS_real.G_k_w(i, i, j, l));

  MOMS_real.A0_w = 0;
  for (int l = 0; l < w_REAL::dmn_size(); l++)
    for (int j = 0; j < k_DCA::dmn_size(); j++)
      for (int i = 0; i < 2 * b::dmn_size(); i++)
        MOMS_real.A0_w(l) -= imag(MOMS_real.G0_k_w(i, i, j, l));

  MOMS_real.A_w *= 1. / double(M_PI * k_DCA::dmn_size() * 2 * b::dmn_size());
  MOMS_real.A0_w *= 1. / double(M_PI * k_DCA::dmn_size() * 2 * b::dmn_size());

  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n" << dca::util::print_time() << std::endl;
  }

  sp_Greens_function_obj.compute_S_k_w(MOMS_imag.G_k_w, MOMS_imag.G0_k_w, MOMS_imag.Sigma);
  sp_Greens_function_obj.compute_S_k_w(MOMS_real.G_k_w, MOMS_real.G0_k_w, MOMS_real.Sigma);

  if (concurrency.id() == concurrency.first()) {
    std::cout << dca::util::print_time() << std::endl;
  }
}

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
template <typename dca_info_struct_t>
void cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::finalize(
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

template <LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
void cluster_solver<ADVANCED_ED_CLUSTER_SOLVER, device_t, parameters_type, MOMS_type>::write(
    std::string file_name) {
  IO::FORMAT FORMAT = parameters.get_output_format();

  std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

  switch (FORMAT) {
    case IO::JSON: {
      IO::writer<IO::JSON> writer;
      {
        writer.open_file(file_name);

        parameters.write(writer);
        MOMS_imag.write(writer);
        MOMS_real.write(writer);

        if (parameters.get_vertex_measurement_type() != NONE) {
          std::cout << "\n\n\t\t start writing tp-Greens-function\n\n";
          tp_Greens_function_obj.write(writer);
        }

        writer.close_file();
      }
    } break;

    case IO::HDF5: {
      IO::writer<IO::HDF5> writer;
      {
        writer.open_file(file_name);

        std::cout << "\n\n\t\t start writing parameters\n\n";
        parameters.write(writer);

        std::cout << "\n\n\t\t start writing MOMS_imag\n\n";
        MOMS_imag.write(writer);

        std::cout << "\n\n\t\t start writing MOMS_real\n\n";
        MOMS_real.write(writer);

        if (parameters.get_vertex_measurement_type() != NONE) {
          std::cout << "\n\n\t\t start writing tp-Greens-function\n\n";
          tp_Greens_function_obj.write(writer);
        }

        writer.close_file();
      }
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_CLUSTER_SOLVER_H
