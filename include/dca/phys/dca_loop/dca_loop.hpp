// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This class executes the DCA(+) loop.

#ifndef DCA_PHYS_DCA_LOOP_DCA_LOOP_HPP
#define DCA_PHYS_DCA_LOOP_DCA_LOOP_HPP

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <variant>

#include "dca/distribution/dist_types.hpp"
#include "dca/function/domains.hpp"
#include "dca/io/filesystem.hpp"
#include "dca/io/writer.hpp"
#include "dca/io/io_types.hpp"
#include "dca/phys/dca_algorithms/compute_greens_function.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_mapping/cluster_exclusion.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_sp.hpp"
#include "dca/phys/dca_step/cluster_mapping/double_counting_correction.hpp"
#include "dca/phys/dca_step/cluster_mapping/update_chemical_potential.hpp"
#include "dca/phys/dca_step/cluster_solver/high_temperature_series_expansion/high_temperature_series_expansion_solver.hpp"
#include "dca/phys/dca_step/lattice_mapping/lattice_mapping_sp.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/util/print_time.hpp"
#include "dca/util/signal_handler.hpp"

namespace dca {
namespace phys {
// dca::phys::

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType,
          DistType DIST = DistType::NONE>
class DcaLoop {
public:
  static constexpr DistType DT = DIST;
  using profiler_type = typename ParametersType::profiler_type;
  using concurrency_type = typename ParametersType::concurrency_type;

  using Lattice = typename ParametersType::lattice_type;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using k_DCA =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::CLUSTER,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  using k_HOST =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::LATTICE_SP,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;

  using cluster_exclusion_type = clustermapping::cluster_exclusion<ParametersType, DcaDataType>;
  using double_counting_correction_type =
      clustermapping::double_counting_correction<ParametersType, DcaDataType>;
  using coarsegraining_sp_type = clustermapping::CoarsegrainingSp<ParametersType>;
  using lattice_map_sp_type = latticemapping::lattice_mapping_sp<ParametersType, k_DCA, k_HOST>;
  using update_chemical_potential_type =
      clustermapping::update_chemical_potential<ParametersType, DcaDataType, coarsegraining_sp_type>;
  using HTS_solver_type =
      solver::HighTemperatureSeriesExpansionSolver<dca::linalg::CPU, ParametersType, DcaDataType>;

  DcaLoop(ParametersType& parameters_ref, DcaDataType& MOMS_ref, concurrency_type& concurrency_ref);

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
  void adjust_impurity_self_energy();

  void perform_cluster_exclusion_step();

  double solve_cluster_problem(int DCA_iteration);

  void perform_lattice_mapping();

  void update_DCA_loop_data_functions(int DCA_iteration);

  void logSelfEnergy(int i);

  ParametersType& parameters;
  DcaDataType& MOMS;
  concurrency_type& concurrency;

private:
  DcaLoopData<ParametersType> DCA_info_struct;

  cluster_exclusion_type cluster_exclusion_obj;
  double_counting_correction_type double_counting_correction_obj;

  coarsegraining_sp_type cluster_mapping_obj;
  lattice_map_sp_type lattice_mapping_obj;

  update_chemical_potential_type update_chemical_potential_obj;

  std::string file_name_;
  std::shared_ptr<io::Writer<concurrency_type>> output_file_;

  unsigned dca_iteration_ = 0;

protected:
  MCIntegratorType monte_carlo_integrator_;
};

/** setup objects with loop scope lifetime
 *
 *  \todo The way output_file_ is handled needs to be changed to support
 *  generalized parellel io
 */
template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::DcaLoop(
    ParametersType& parameters_ref, DcaDataType& MOMS_ref, concurrency_type& concurrency_ref)
    : parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(concurrency_ref),
      cluster_exclusion_obj(parameters, MOMS),
      double_counting_correction_obj(parameters, MOMS),
      cluster_mapping_obj(parameters),
      lattice_mapping_obj(parameters),
      update_chemical_potential_obj(parameters, MOMS, cluster_mapping_obj),
#ifdef DCA_HAVE_ADIOS2
      output_file_(
          parameters.get_output_format() == "ADIOS2"
              ? std::make_shared<io::Writer<concurrency_type>>(
                    concurrency.get_adios(), concurrency_ref, parameters.get_output_format(), false)
              : std::make_shared<io::Writer<concurrency_type>>(
                           concurrency.get_adios(), concurrency_ref, parameters.get_output_format(),
                           false)),
#else
      output_file_(std::make_shared<io::Writer<concurrency_type>>(
								  concurrency_ref, parameters.get_output_format(), false)),
#endif
      monte_carlo_integrator_(parameters_ref, MOMS_ref, output_file_) {
  file_name_ = parameters.get_directory() + parameters.get_filename_dca();

  if (concurrency.id() == concurrency.first()) {
    // dca::util::SignalHandler<concurrency_type>::registerFile(output_file_);
    std::cout << "\n\n\t" << __FUNCTION__ << " has started \t" << dca::util::print_time() << "\n\n";
  }
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::write() {
  // We assume DCALoop write is called once at the end of the run and we leave a step open for it to
  // write into.

  //output_file_->begin_step();

  if (concurrency.id() == concurrency.first()) {
    // This should probably happen first not at the end
    parameters.write(*output_file_);
    MOMS.write(*output_file_);
    monte_carlo_integrator_.write(*output_file_);
    DCA_info_struct.write(*output_file_, concurrency);
    // None of this kludgy file shuffling with ADIOS2 each iteration is a step in the ADIOS2 output
  }
  output_file_->end_step();

  if (concurrency.id() == concurrency.first()) {
    if (output_file_ && !(output_file_->isADIOS2())) {
      output_file_->close_file();
      output_file_.reset();
      std::error_code code;
      filesystem::rename(file_name_ + ".tmp", file_name_, code);
      if (code) {
        std::cerr << "Failed to rename file." << std::endl;
      }
    }
  }
  // For ADIOS2 every rank has an output_file_ so close them
  // data.
  if (output_file_ != nullptr && output_file_->isADIOS2()) {
    // if (concurrency.id() == concurrency.first()) {
    output_file_->close_file();
    std::cout << "ADIOS2 output closed\n";
    //}
  }
}

template <typename ParametersType, typename DDT, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DDT, MCIntegratorType, DIST>::initialize() {
  static_assert(std::is_same<DDT, dca::phys::DcaData<ParametersType, DIST>>::value);
  int last_completed = -1;
  auto& autoresume_filename = parameters.get_autoresume_filename();
  io::IOType iotype = io::extensionToIOType(autoresume_filename);
  if (parameters.autoresume()) {
#ifdef DCA_HAVE_ADIOS2
    if (iotype == io::IOType::ADIOS2)
      last_completed = DCA_info_struct.readData(autoresume_filename, parameters.get_output_format(),
                                                concurrency, concurrency.get_adios());
    else
#endif
      last_completed =
          DCA_info_struct.readData(autoresume_filename, parameters.get_output_format(), concurrency);

    if (last_completed >= 0) {
      if (concurrency.id() == concurrency.first())
        std::cout << "\n   *******  Resuming DCA from iteration " << last_completed + 1
                  << "  *******\n"
                  << std::endl;

      dca_iteration_ = std::min(last_completed + 1, parameters.get_dca_iterations() - 1);
#ifdef DCA_HAVE_ADIOS2
      if (iotype == io::IOType::ADIOS2)
        MOMS.initializeSigma(concurrency.get_adios(), autoresume_filename);
      else
#endif
        MOMS.initializeSigma(autoresume_filename);
      perform_lattice_mapping();
    }
  }
  else if (parameters.get_initial_self_energy() != "zero") {
#ifdef DCA_HAVE_ADIOS2
    if (io::extensionToIOType(parameters.get_initial_self_energy()) == io::IOType::ADIOS2)
      MOMS.initializeSigma(concurrency.get_adios(), parameters.get_initial_self_energy());
    else
#endif
      MOMS.initializeSigma(parameters.get_initial_self_energy());

    perform_lattice_mapping();
  }

  // In ADIOS2 every rank has an output_file_ for the purposes of writing per rank
  // or per configuration data.
  if (output_file_ && output_file_->isADIOS2()) {
    output_file_->open_file(file_name_, parameters.autoresume() ? false : true);
  }
  else {
    if (concurrency.id() == concurrency.first())
      output_file_->open_file(file_name_ + ".tmp", parameters.autoresume() ? false : true);
  }
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::execute() {
  for (; dca_iteration_ < parameters.get_dca_iterations(); dca_iteration_++) {
    output_file_->begin_step();

    adjust_chemical_potential();

    perform_cluster_mapping();

    adjust_coarsegrained_self_energy();  // double-counting-correction

    perform_cluster_exclusion_step();

    double L2_Sigma_difference =
        solve_cluster_problem(dca_iteration_);  // returned from cluster_solver::finalize

    adjust_impurity_self_energy();  // double-counting-correction

    perform_lattice_mapping();

    update_DCA_loop_data_functions(dca_iteration_);  // Really just updates

    // This writes the self energy and DCA_loop_data.  Sort of a check point.
    logSelfEnergy(dca_iteration_);

    if (L2_Sigma_difference <
        parameters.get_dca_accuracy()) {  // set the acquired accuracy on |Sigma_QMC - Sigma_cg|
      break;
    }

    // As long as this isn't the last iteration where this is handled by the finalize we want these here
    if (dca_iteration_ != parameters.get_dca_iterations() - 1) {
      if (output_file_ && (output_file_->isADIOS2() || output_file_->isHDF5())) {
        if (parameters.dump_every_iteration()) {
          if (concurrency.id() == concurrency.first()) {
            // This normally gets done in finalize before the dump.
            perform_cluster_mapping_self_energy();
            MOMS.compute_Sigma_bands();
            MOMS.compute_single_particle_properties();
            MOMS.write(*output_file_);
          }
        }
      }
      output_file_->end_step();
    }

    if (parameters.do_not_update_sigma()) {
	if (parameters.get_initial_self_energy() == "zero")
	  throw std::runtime_error("If there is no initial self energy it must be updated!");
#ifdef DCA_HAVE_ADIOS2
	if (io::extensionToIOType(parameters.get_initial_self_energy()) == io::IOType::ADIOS2)
	  MOMS.initializeSigma(concurrency.get_adios(), parameters.get_initial_self_energy());
	else
#endif
	  MOMS.initializeSigma(parameters.get_initial_self_energy());
    }

  }
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::finalize() {
  perform_cluster_mapping_self_energy();
  MOMS.compute_Sigma_bands();
  MOMS.compute_single_particle_properties();
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::adjust_chemical_potential() {
  if (parameters.adjust_chemical_potential())
    update_chemical_potential_obj.execute();
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::perform_cluster_mapping() {
  perform_cluster_mapping_self_energy();
  perform_cluster_mapping_Greens_function();
  // perform_cluster_exclusion_step();
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::perform_cluster_mapping_self_energy() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t coarsegrain-Selfenergy " << dca::util::print_time();

  profiler_type profiler("coarsegrain-Selfenergy", "DCA", __LINE__);

  if (parameters.do_dca_plus())
    cluster_mapping_obj.compute_S_K_w(MOMS.Sigma_lattice, MOMS.Sigma_cluster);
  else
    MOMS.Sigma_cluster = MOMS.Sigma;

  MOMS.print_Sigma_QMC_versus_Sigma_cg();

  Symmetrize<ParametersType>::execute(MOMS.Sigma_cluster, MOMS.H_symmetry);
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType,
             DIST>::perform_cluster_mapping_Greens_function() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t coarsegrain-Greens-function " << dca::util::print_time();

  profiler_type profiler("coarsegrain-Greens-function", "DCA", __LINE__);

  // Finite-size QMC
  if (parameters.do_finite_size_qmc())
    compute_G_k_w(MOMS.H_DCA, MOMS.Sigma, parameters.get_chemical_potential(),
                  parameters.get_coarsegraining_threads(), MOMS.G_k_w);
  // DCA+
  else if (parameters.do_dca_plus())
    cluster_mapping_obj.compute_G_K_w(MOMS.Sigma_lattice, MOMS.G_k_w);
  // Standard DCA
  else
    cluster_mapping_obj.compute_G_K_w(MOMS.Sigma, MOMS.G_k_w);

  Symmetrize<ParametersType>::execute(MOMS.G_k_w, MOMS.H_symmetry);
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::adjust_coarsegrained_self_energy() {
  double_counting_correction_obj.execute_before_solver();
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::adjust_impurity_self_energy() {
  double_counting_correction_obj.execute_after_solver();
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::perform_cluster_exclusion_step() {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t cluster-exclusion-step " << dca::util::print_time();

  profiler_type profiler("cluster-exclusion-step", "DCA", __LINE__);

  cluster_exclusion_obj.execute();
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
double DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::solve_cluster_problem(
    int DCA_iteration) {
  //  static_assert(std::is_same<DcaDataType, ::DcaDataType<DIST>>::value);
  // static_assert(std::is_same<MCIntegratorType, ::ClusterSolver<DIST>>::value);
  {
    profiler_type profiler("initialize cluster-solver", "DCA", __LINE__);
    monte_carlo_integrator_.initialize(DCA_iteration);
  }

  {
    profiler_type profiler("Quantum Monte Carlo integration", "DCA", __LINE__);
    monte_carlo_integrator_.integrate();
  }

  if (output_file_ && output_file_->isADIOS2())
    output_file_->flush();

  {
    if (concurrency.id() == concurrency.first())
      std::cout << "start Monte Carlo integration finalize.\n";

    profiler_type profiler("finalize cluster-solver", "DCA", __LINE__);
    double L2_Sigma_difference = monte_carlo_integrator_.finalize(DCA_info_struct);
    if (concurrency.id() == concurrency.first())
      std::cout << "Monte Carlo integration finalized.\n";

    return L2_Sigma_difference;
  }
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::perform_lattice_mapping() {
  profiler_type profiler("lattice-mapping", "DCA", __LINE__);

  if (concurrency.id() == concurrency.first())
    std::cout << "\n\t\t lattice-mapping " << dca::util::print_time();

  if (parameters.do_dca_plus() || parameters.doPostInterpolation()) {
    if (parameters.hts_approximation()) {
      DcaDataType MOMS_HTS(parameters);

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

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::update_DCA_loop_data_functions(int i) {
  DCA_info_struct.density(i) = update_chemical_potential_obj.compute_density();
  DCA_info_struct.chemical_potential(i) = parameters.get_chemical_potential();

  if (concurrency.id() == concurrency.first()) {
    std::cout << "\n\n\t\t\t total-density : " << DCA_info_struct.density(i)
              << "\t (time : " << dca::util::print_time() << ")\n\n";
  }

  for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++)
    DCA_info_struct.orbital_occupancies(l1, i) = MOMS.orbital_occupancy(l1);

  for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++)
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
      DCA_info_struct.n_k(l1, k_ind, i) = 1. - std::abs(MOMS.G_k_t(l1, l1, k_ind, 0));

  for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++)
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
      // TODO: Use t::dmn_size() instead of parameters.get_sp_time_intervals().
      DCA_info_struct.A_k(l1, k_ind, i) =
          std::abs(MOMS.G_k_t(l1, l1, k_ind, parameters.get_sp_time_intervals() / 2)) *
          parameters.get_beta() / M_PI;
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::logSelfEnergy(int i) {
  DCA_info_struct.last_completed_iteration = i;

  if (output_file_ && concurrency.id() == concurrency.first()) {
    output_file_->open_group("functions");
    output_file_->execute(MOMS.Sigma);
    output_file_->close_group();

    output_file_->open_group("parameters");
    output_file_->open_group("physics");
    output_file_->execute("chemical-potential", parameters.get_chemical_potential());
    output_file_->close_group();
    output_file_->close_group();

    DCA_info_struct.write(*output_file_, concurrency);
  }
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_LOOP_DCA_LOOP_HPP
