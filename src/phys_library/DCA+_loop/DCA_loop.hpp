// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//
// This class executes the DCA(+) loop.

#ifndef PHYS_LIBRARY_DCA_LOOP_DCA_LOOP_HPP
#define PHYS_LIBRARY_DCA_LOOP_DCA_LOOP_HPP

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "dca/phys/DCA_step/cluster_mapping/cluster_exclusion.hpp"
#include "dca/phys/DCA_step/cluster_mapping/double_counting_correction.hpp"
#include "dca/phys/DCA_step/cluster_mapping/update_chemical_potential.hpp"
#include "dca/util/print_time.hpp"
#include "comp_library/function_library/include_function_library.h"
#include "comp_library/IO_library/IO.hpp"
#include "phys_library/DCA+_loop/DCA_loop_data.hpp"
#include "phys_library/DCA+_step/cluster_mapping/coarsegraining_step/coarsegraining_sp.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/high_temperature_series_expansion_solver.h"
#include "phys_library/DCA+_step/lattice_mapping/lattice_mapping_sp.h"
#include "phys_library/DCA+_step/symmetrization/symmetrize.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"

namespace DCA {

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
class DCA_loop {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using k_DCA = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, CLUSTER,
                                     MOMENTUM_SPACE, BRILLOUIN_ZONE>>;
  using k_HOST = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, LATTICE_SP,
                                      MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  using cluster_exclusion_type = DCA::cluster_exclusion<parameters_type, MOMS_type>;
  using double_counting_correction_type = DCA::double_counting_correction<parameters_type, MOMS_type>;
  using coarsegraining_sp_type = DCA::coarsegraining_sp<parameters_type, k_DCA>;
  using lattice_map_sp_type = DCA::lattice_mapping_sp<parameters_type, k_DCA, k_HOST>;
  using update_chemical_potential_type =
      DCA::update_chemical_potential<parameters_type, MOMS_type, coarsegraining_sp_type>;
  using HTS_solver_type =
      DCA::cluster_solver<DCA::HIGH_TEMPERATURE_SERIES, dca::linalg::CPU, parameters_type, MOMS_type>;

public:
  DCA_loop(parameters_type& parameters_ref, MOMS_type& MOMS_ref, concurrency_type& concurrency_ref);

  ~DCA_loop();

  void read();

  void write();

  void initialize();

  void execute();

  void finalize();

protected:
  void adjust_chemical_potential();

  void perform_cluster_mapping();

  void perform_cluster_mapping_self_energy();

  void perform_cluster_mapping_Greens_function();

  void adjust_coarsegrained_self_energy();

  void perform_cluster_exclusion_step();

  double solve_cluster_problem(int DCA_iteration);

  void adjust_impurity_self_energy();

  void perform_lattice_mapping();

  void perform_lattice_mapping_with_HTS();
  void perform_lattice_mapping_without_HTS();

  void update_DCA_loop_data_functions(int DCA_iteration);

protected:
  parameters_type& parameters;
  MOMS_type& MOMS;
  concurrency_type& concurrency;

private:
  DCA_loop_data<parameters_type> DCA_info_struct;

  cluster_exclusion_type cluster_exclusion_obj;
  double_counting_correction_type double_counting_correction_obj;

  coarsegraining_sp_type cluster_mapping_obj;
  lattice_map_sp_type lattice_mapping_obj;

  update_chemical_potential_type update_chemical_potential_obj;

protected:
  Monte_Carlo_Integrator_type monte_carlo_integrator_;
};

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::DCA_loop(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref, concurrency_type& concurrency_ref)
    : parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(concurrency_ref),

      DCA_info_struct(),

      cluster_exclusion_obj(parameters, MOMS),
      double_counting_correction_obj(parameters, MOMS),

      cluster_mapping_obj(parameters),
      lattice_mapping_obj(parameters),

      update_chemical_potential_obj(parameters, MOMS, cluster_mapping_obj),

      monte_carlo_integrator_(parameters_ref, MOMS_ref) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t" << __FUNCTION__ << " has started \t" << dca::util::print_time() << "\n\n";
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::~DCA_loop() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t" << __FUNCTION__ << " has finished \t" << dca::util::print_time()
              << "\n\n";
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::read() {
  if (parameters.get_Sigma_file() != "zero")
    MOMS.read(parameters.get_Sigma_file());
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::write() {
  IO::FORMAT FORMAT = parameters.get_output_format();
  std::string file_name = parameters.get_directory() + parameters.get_output_file_name();

  std::cout << "\n\n\t\t start writing " << file_name << "\t" << dca::util::print_time() << "\n\n";

  switch (FORMAT) {
    case IO::JSON: {
      IO::writer<IO::JSON> writer;
      {
        writer.open_file(file_name);

        parameters.write(writer);
        MOMS.write(writer);
        monte_carlo_integrator_.write(writer);
        DCA_info_struct.write(writer);

        writer.close_file();
      }
    } break;

    case IO::HDF5: {
      IO::writer<IO::HDF5> writer;
      {
        writer.open_file(file_name);

        parameters.write(writer);
        MOMS.write(writer);
        monte_carlo_integrator_.write(writer);
        DCA_info_struct.write(writer);

        writer.close_file();
      }
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::initialize() {
  if (parameters.get_Sigma_file() != "zero") {
    MOMS.initialize_Sigma();

    perform_lattice_mapping();
  }
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::execute() {
  for (int i = 0; i < parameters.get_DCA_iterations(); i++) {
    adjust_chemical_potential();

    perform_cluster_mapping();

    adjust_coarsegrained_self_energy();  // double-counting-correction

    perform_cluster_exclusion_step();

    double L2_Sigma_difference =
        solve_cluster_problem(i);  // returned from cluster_solver::finalize

    adjust_impurity_self_energy();  // double-counting-correction

    perform_lattice_mapping();

    update_DCA_loop_data_functions(i);

    if (L2_Sigma_difference <
        parameters.get_DCA_accuracy())  // set the acquired accuracy on |Sigma_QMC - Sigma_cg|
      break;
  }
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::finalize() {
  perform_cluster_mapping_self_energy();

  MOMS.compute_Sigma_bands();

  MOMS.compute_single_particle_properties();
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::adjust_chemical_potential() {
  if (parameters.adjust_chemical_potential())
    update_chemical_potential_obj.execute();
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::perform_cluster_mapping() {
  perform_cluster_mapping_self_energy();

  perform_cluster_mapping_Greens_function();

  // perform_cluster_exclusion_step();
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type,
              Monte_Carlo_Integrator_type>::perform_cluster_mapping_self_energy() {
  if (concurrency.id() == 0)
    std::cout << "\n\t\t coarsegrain-Selfenergy " << dca::util::print_time();

  profiler_type profiler("coarsegrain-Selfenergy", "DCA", __LINE__);

  if (parameters.do_DCA_plus())
    cluster_mapping_obj.compute_S_K_w(MOMS.Sigma_lattice, MOMS.Sigma_cluster);
  else
    MOMS.Sigma_cluster = MOMS.Sigma;

  MOMS.print_Sigma_QMC_versus_Sigma_cg();

  symmetrize::execute(MOMS.Sigma_cluster, MOMS.H_symmetry);
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type,
              Monte_Carlo_Integrator_type>::perform_cluster_mapping_Greens_function() {
  if (concurrency.id() == 0)
    std::cout << "\n\t\t coarsegrain-Greens-function " << dca::util::print_time();

  profiler_type profiler("coarsegrain-Greens-function", "DCA", __LINE__);

  if (parameters.do_DCA_plus())
    cluster_mapping_obj.compute_G_K_w(MOMS.H_HOST, MOMS.Sigma_lattice, MOMS.G_k_w);
  else
    cluster_mapping_obj.compute_G_K_w(MOMS.H_HOST, MOMS.Sigma, MOMS.G_k_w);

  symmetrize::execute(MOMS.G_k_w, MOMS.H_symmetry);
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::adjust_coarsegrained_self_energy() {
  double_counting_correction_obj.execute_before_solver();
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::perform_cluster_exclusion_step() {
  if (concurrency.id() == 0)
    std::cout << "\n\t\t cluster-exclusion-step " << dca::util::print_time();

  profiler_type profiler("cluster-exclusion-step", "DCA", __LINE__);

  cluster_exclusion_obj.execute();
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
double DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::solve_cluster_problem(
    int DCA_iteration) {
  {
    profiler_type profiler("initialize cluster-solver", "DCA", __LINE__);
    monte_carlo_integrator_.initialize(DCA_iteration);
  }

  {
    profiler_type profiler("Quantum Monte Carlo integration", "DCA", __LINE__);
    monte_carlo_integrator_.integrate();
  }

  {
    profiler_type profiler("finalize cluster-solver", "DCA", __LINE__);
    double L2_Sigma_difference = monte_carlo_integrator_.finalize(DCA_info_struct);

    return L2_Sigma_difference;
  }
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::adjust_impurity_self_energy() {
  double_counting_correction_obj.execute_after_solver();
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::perform_lattice_mapping() {
  profiler_type profiler("lattice-mapping", "DCA", __LINE__);

  if (concurrency.id() == 0)
    std::cout << "\n\t\t lattice-mapping " << dca::util::print_time();

  if (parameters.do_DCA_plus()) {
    if (parameters.use_HTS_approximation()) {
      MOMS_type MOMS_HTS(parameters);

      MOMS_HTS.H_HOST = MOMS.H_HOST;
      MOMS_HTS.H_interactions = MOMS.H_interactions;

      HTS_solver_type HTS_solver(parameters, MOMS_HTS);

      lattice_mapping_obj.execute_with_HTS_approximation(
          MOMS_HTS, HTS_solver, cluster_mapping_obj, MOMS.Sigma, MOMS.Sigma_lattice_interpolated,
          MOMS.Sigma_lattice_coarsegrained, MOMS.Sigma_lattice);
    }
    else {
      lattice_mapping_obj.execute(MOMS.Sigma, MOMS.Sigma_lattice_interpolated,
                                  MOMS.Sigma_lattice_coarsegrained, MOMS.Sigma_lattice);
    }
  }
}

template <class parameters_type, class MOMS_type, class Monte_Carlo_Integrator_type>
void DCA_loop<parameters_type, MOMS_type, Monte_Carlo_Integrator_type>::update_DCA_loop_data_functions(
    int i) {
  DCA_info_struct.density(i) = update_chemical_potential_obj.compute_density();
  DCA_info_struct.chemical_potential(i) = parameters.get_chemical_potential();

  if (concurrency.id() == 0) {
    std::cout << "\n\n\t\t\t total-density : " << DCA_info_struct.density(i)
              << "\t (time : " << dca::util::print_time() << ")\n\n";
  }

  for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++)
    DCA_info_struct.orbital_occupancies(l1, i) = MOMS.orbital_occupancy(l1);

  for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++)
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
      DCA_info_struct.n_k(l1, k_ind, i) = 1. - MOMS.G_k_t(l1, l1, k_ind, 0);

  for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++)
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
      DCA_info_struct.A_k(l1, k_ind, i) =
          MOMS.G_k_t(l1, l1, k_ind, parameters.get_sp_time_intervals() / 2) *
          parameters.get_beta() / M_PI;
}
}

#endif  // PHYS_LIBRARY_DCA_LOOP_DCA_LOOP_HPP
