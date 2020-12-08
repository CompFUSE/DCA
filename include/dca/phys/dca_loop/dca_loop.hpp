// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//
// This class executes the DCA(+) loop.

#ifndef DCA_PHYS_DCA_LOOP_DCA_LOOP_HPP
#define DCA_PHYS_DCA_LOOP_DCA_LOOP_HPP

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "dca/distribution/dist_types.hpp"
#include "dca/function/domains.hpp"
#include "dca/io/filesystem.hpp"
#include "dca/io/writer.hpp"
#include "dca/phys/dca_algorithms/compute_greens_function.hpp"
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

  DcaLoop(ParametersType& parameters_ref, DcaDataType& data__ref, concurrency_type& concurrency_ref);

  void write();

  void initialize();
  void execute();
  void finalize();

private:
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

  void readInitialStatus(const std::string& filename);

  ParametersType& parameters_;
  DcaDataType& data_;
  concurrency_type& concurrency_;

  DcaLoopData<ParametersType> DCA_info_struct;

  cluster_exclusion_type cluster_exclusion_obj;
  double_counting_correction_type double_counting_correction_obj;

  coarsegraining_sp_type cluster_mapping_obj;
  lattice_map_sp_type lattice_mapping_obj;

  update_chemical_potential_type update_chemical_potential_obj;

  std::string file_name_;
  std::shared_ptr<io::Writer> output_file_;

  unsigned dca_iteration_ = 0;

protected:
  MCIntegratorType monte_carlo_integrator_;
};

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::DcaLoop(
    ParametersType& parameters_ref, DcaDataType& data_ref, concurrency_type& concurrency_ref)
    : parameters_(parameters_ref),
      data_(data_ref),
      concurrency_(concurrency_ref),

      DCA_info_struct(),

      cluster_exclusion_obj(parameters_, data_),
      double_counting_correction_obj(parameters_, data_),

      cluster_mapping_obj(parameters_),
      lattice_mapping_obj(parameters_),

      update_chemical_potential_obj(parameters_, data_, cluster_mapping_obj),

      output_file_(concurrency_.id() == concurrency_.first()
                       ? std::make_shared<io::Writer>(parameters_.get_output_format(), false)
                       : nullptr),

      monte_carlo_integrator_(parameters_ref, data_, output_file_) {
  dca::util::SignalHandler::registerFile(output_file_);

  std::cout << "\n\n\t" << __FUNCTION__ << " has started \t" << dca::util::print_time() << "\n\n";
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::write() {
  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "\n\n\t\t start writing " << file_name_ << "\t" << dca::util::print_time() << "\n\n";

    output_file_->set_verbose(true);

    parameters_.write(*output_file_);
    data_.write(*output_file_);
    monte_carlo_integrator_.write(*output_file_);
    DCA_info_struct.write(*output_file_);

    output_file_->close_file();
    output_file_.reset();

    std::error_code code;
    filesystem::rename(file_name_ + ".tmp", file_name_, code);
    if (code) {
      std::cerr << "Failed to rename file." << std::endl;
    }
  }
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::initialize() {
  int last_completed = -1;

  if (parameters_.autoresume()) {  // Try to read state of previous run.
    last_completed = DCA_info_struct.tryToRead(file_name_ + ".tmp", parameters_.get_output_format(),
                                               concurrency_);
  }
  if (last_completed >= 0) {
    if (concurrency_.id() == concurrency_.first())
      std::cout << "\n   *******  Resuming DCA from iteration " << last_completed + 1 << "  *******\n"
                << std::endl;

    dca_iteration_ = std::min(last_completed + 1, parameters_.get_dca_iterations() - 1);
    readInitialStatus(file_name_ + ".tmp");
  }
  else if (parameters_.get_initial_self_energy() != "zero") {
    readInitialStatus(parameters_.get_initial_self_energy());
  }

  if (concurrency_.id() == concurrency_.first()) {
    output_file_->open_file(file_name_ + ".tmp", parameters_.autoresume() ? false : true);
  }
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::readInitialStatus(
    const std::string& filename) {
  io::Buffer buffer;
  bool buffer_read = false;

  if (concurrency_.id() == concurrency_.first()) {
    io::HDF5Reader reader;
    reader.open_file(filename);

    if (parameters_.adjust_chemical_potential()) {
      reader.open_group("parameters");
      reader.open_group("physics");
      reader.execute("chemical-potential", parameters_.get_chemical_potential());
      reader.close_group();
      reader.close_group();
    }

    reader.open_group("functions");
    reader.execute(data_.Sigma);
    reader.close_group();

    if (parameters_.store_configuration()) {
      reader.open_group("Configurations");
      buffer_read = reader.execute("sample", buffer);
      reader.close_group();
    }
  }

  concurrency_.broadcast(parameters_.get_chemical_potential());
  concurrency_.broadcast(data_.Sigma);

  concurrency_.broadcast(buffer_read);
  if (buffer_read) {
    concurrency_.broadcast(buffer);
    monte_carlo_integrator_.setSampleConfiguration(buffer);
  }

  perform_lattice_mapping();
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::execute() {
  for (; dca_iteration_ < parameters_.get_dca_iterations(); dca_iteration_++) {
    adjust_chemical_potential();

    perform_cluster_mapping();

    adjust_coarsegrained_self_energy();  // double-counting-correction

    perform_cluster_exclusion_step();

    double L2_Sigma_difference =
        solve_cluster_problem(dca_iteration_);  // returned from cluster_solver::finalize

    adjust_impurity_self_energy();  // double-counting-correction

    perform_lattice_mapping();

    update_DCA_loop_data_functions(dca_iteration_);
    logSelfEnergy(dca_iteration_);  // Write a check point.

    if (L2_Sigma_difference <
        parameters_.get_dca_accuracy())  // set the acquired accuracy on |Sigma_QMC - Sigma_cg|
      break;
  }
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::finalize() {
  perform_cluster_mapping_self_energy();
  data_.compute_Sigma_bands();
  data_.compute_single_particle_properties();
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::adjust_chemical_potential() {
  if (parameters_.adjust_chemical_potential())
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
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t coarsegrain-Selfenergy " << dca::util::print_time();

  profiler_type profiler("coarsegrain-Selfenergy", "DCA", __LINE__);

  if (parameters_.do_dca_plus())
    cluster_mapping_obj.compute_S_K_w(data_.Sigma_lattice, data_.Sigma_cluster);
  else
    data_.Sigma_cluster = data_.Sigma;

  data_.print_Sigma_QMC_versus_Sigma_cg();

  symmetrize::execute<Lattice>(data_.Sigma_cluster, data_.H_symmetry);
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType,
             DIST>::perform_cluster_mapping_Greens_function() {
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t coarsegrain-Greens-function " << dca::util::print_time();

  profiler_type profiler("coarsegrain-Greens-function", "DCA", __LINE__);

  // Finite-size QMC
  if (parameters_.do_finite_size_qmc())
    compute_G_k_w(data_.H_DCA, data_.Sigma, parameters_.get_chemical_potential(),
                  parameters_.get_coarsegraining_threads(), data_.G_k_w);
  // DCA+
  else if (parameters_.do_dca_plus())
    cluster_mapping_obj.compute_G_K_w(data_.Sigma_lattice, data_.G_k_w);
  // Standard DCA
  else
    cluster_mapping_obj.compute_G_K_w(data_.Sigma, data_.G_k_w);

  symmetrize::execute<Lattice>(data_.G_k_w, data_.H_symmetry);
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
  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t cluster-exclusion-step " << dca::util::print_time();

  profiler_type profiler("cluster-exclusion-step", "DCA", __LINE__);

  cluster_exclusion_obj.execute();
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
double DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::solve_cluster_problem(
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

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::perform_lattice_mapping() {
  profiler_type profiler("lattice-mapping", "DCA", __LINE__);

  if (concurrency_.id() == concurrency_.first())
    std::cout << "\n\t\t lattice-mapping " << dca::util::print_time();

  if (parameters_.do_dca_plus() || parameters_.doPostInterpolation()) {
    if (parameters_.hts_approximation()) {
      DcaDataType data__HTS(parameters_);

      data__HTS.H_HOST = data_.H_HOST;
      data__HTS.H_interactions = data_.H_interactions;

      HTS_solver_type HTS_solver(parameters_, data__HTS);

      lattice_mapping_obj.execute_with_HTS_approximation(
          data__HTS, HTS_solver, cluster_mapping_obj, data_.Sigma, data_.Sigma_lattice_interpolated,
          data_.Sigma_lattice_coarsegrained, data_.Sigma_lattice);
    }
    else {
      lattice_mapping_obj.execute(data_.Sigma, data_.Sigma_lattice_interpolated,
                                  data_.Sigma_lattice_coarsegrained, data_.Sigma_lattice);
    }
  }
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::update_DCA_loop_data_functions(int i) {
  DCA_info_struct.density(i) = update_chemical_potential_obj.compute_density();
  DCA_info_struct.chemical_potential(i) = parameters_.get_chemical_potential();

  if (concurrency_.id() == concurrency_.first()) {
    std::cout << "\n\n\t\t\t total-density : " << DCA_info_struct.density(i)
              << "\t (time : " << dca::util::print_time() << ")\n\n";
  }

  for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++)
    DCA_info_struct.orbital_occupancies(l1, i) = data_.orbital_occupancy(l1);

  for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++)
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
      DCA_info_struct.n_k(l1, k_ind, i) = 1. - std::abs(data_.G_k_t(l1, l1, k_ind, 0));

  for (int l1 = 0; l1 < b::dmn_size() * s::dmn_size(); l1++)
    for (int k_ind = 0; k_ind < k_DCA::dmn_size(); k_ind++)
      // TODO: Use t::dmn_size() instead of parameters_.get_sp_time_intervals().
      DCA_info_struct.A_k(l1, k_ind, i) =
          std::abs(data_.G_k_t(l1, l1, k_ind, parameters_.get_sp_time_intervals() / 2)) *
          parameters_.get_beta() / M_PI;
}

template <typename ParametersType, typename DcaDataType, typename MCIntegratorType, DistType DIST>
void DcaLoop<ParametersType, DcaDataType, MCIntegratorType, DIST>::logSelfEnergy(int i) {
  DCA_info_struct.last_completed_iteration = i;

  if (output_file_) {
    output_file_->open_group("functions");
    output_file_->execute(data_.Sigma);
    output_file_->close_group();

    output_file_->open_group("parameters");
    output_file_->open_group("physics");
    output_file_->execute("chemical-potential", parameters_.get_chemical_potential());
    output_file_->close_group();
    output_file_->close_group();

    DCA_info_struct.write(*output_file_);
  }
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_LOOP_DCA_LOOP_HPP
