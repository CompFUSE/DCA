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
// This class reads, stores, and writes the Monte Carlo solver parameters.
// It is templated on the MC solver name and only implemented when specialized for CT-AUX and
// SS-CT-HYB.

#ifndef DCA_PHYS_PARAMETERS_MC_SOLVER_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_MC_SOLVER_PARAMETERS_HPP

#include <array>
#include <iostream>
#include <stdexcept>
#include "dca/phys/dca_step/cluster_solver/cluster_solver_name.hpp"

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

// Empty class template
template <solver::ClusterSolverName solver_name>
class McSolverParameters {};

// Specialization for CT-AUX
template <>
class McSolverParameters<solver::CT_AUX> {
public:
  McSolverParameters()
      : expansion_parameter_K_(1.),
        initial_configuration_size_(10),
        initial_matrix_size_(128),
        max_submatrix_size_(128),
        neglect_bennett_updates_(false),
        additional_time_measurements_(false) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_expansion_parameter_K() const {
    return expansion_parameter_K_;
  }
  int get_initial_configuration_size() const {
    return initial_configuration_size_;
  }
  int get_initial_matrix_size() const {
    return initial_matrix_size_;
  }
  int get_max_submatrix_size() const {
    return max_submatrix_size_;
  }
  bool neglect_bennett_updates() const {
    return neglect_bennett_updates_;
  }
  bool additional_time_measurements() const {
    return additional_time_measurements_;
  }

private:
  double expansion_parameter_K_;
  int initial_configuration_size_;
  int initial_matrix_size_;
  int max_submatrix_size_;
  bool neglect_bennett_updates_;
  bool additional_time_measurements_;
};

template <typename Concurrency>
int McSolverParameters<solver::CT_AUX>::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(expansion_parameter_K_);
  buffer_size += concurrency.get_buffer_size(initial_configuration_size_);
  buffer_size += concurrency.get_buffer_size(initial_matrix_size_);
  buffer_size += concurrency.get_buffer_size(max_submatrix_size_);
  buffer_size += concurrency.get_buffer_size(neglect_bennett_updates_);
  buffer_size += concurrency.get_buffer_size(additional_time_measurements_);

  return buffer_size;
}

template <typename Concurrency>
void McSolverParameters<solver::CT_AUX>::pack(const Concurrency& concurrency, char* buffer,
                                              int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, expansion_parameter_K_);
  concurrency.pack(buffer, buffer_size, position, initial_configuration_size_);
  concurrency.pack(buffer, buffer_size, position, initial_matrix_size_);
  concurrency.pack(buffer, buffer_size, position, max_submatrix_size_);
  concurrency.pack(buffer, buffer_size, position, neglect_bennett_updates_);
  concurrency.pack(buffer, buffer_size, position, additional_time_measurements_);
}

template <typename Concurrency>
void McSolverParameters<solver::CT_AUX>::unpack(const Concurrency& concurrency, char* buffer,
                                                int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, expansion_parameter_K_);
  concurrency.unpack(buffer, buffer_size, position, initial_configuration_size_);
  concurrency.unpack(buffer, buffer_size, position, initial_matrix_size_);
  concurrency.unpack(buffer, buffer_size, position, max_submatrix_size_);
  concurrency.unpack(buffer, buffer_size, position, neglect_bennett_updates_);
  concurrency.unpack(buffer, buffer_size, position, additional_time_measurements_);
}

// TODO: None of the open_group methods seem to throw.
template <typename ReaderOrWriter>
void McSolverParameters<solver::CT_AUX>::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("CT-AUX");

    try {
      reader_or_writer.execute("expansion-parameter-K", expansion_parameter_K_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("initial-configuration-size", initial_configuration_size_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("initial-matrix-size", initial_matrix_size_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("max-submatrix-size", max_submatrix_size_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("neglect-Bennett-updates", neglect_bennett_updates_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("additional-time-measurements", additional_time_measurements_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}

// Specialization for SS-CT-HYB
template <>
class McSolverParameters<solver::SS_CT_HYB> {
public:
  McSolverParameters()
      : self_energy_tail_cutoff_(0), steps_per_sweep_(0.5), shifts_per_sweep_(0.5) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  int get_self_energy_tail_cutoff() const {
    return self_energy_tail_cutoff_;
  }
  double get_steps_per_sweep() const {
    return steps_per_sweep_;
  }
  double get_shifts_per_sweep() const {
    return shifts_per_sweep_;
  }

private:
  int self_energy_tail_cutoff_;
  double steps_per_sweep_;
  double shifts_per_sweep_;
};

template <typename Concurrency>
int McSolverParameters<solver::SS_CT_HYB>::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(self_energy_tail_cutoff_);
  buffer_size += concurrency.get_buffer_size(steps_per_sweep_);
  buffer_size += concurrency.get_buffer_size(shifts_per_sweep_);

  return buffer_size;
}

template <typename Concurrency>
void McSolverParameters<solver::SS_CT_HYB>::pack(const Concurrency& concurrency, char* buffer,
                                                 int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, self_energy_tail_cutoff_);
  concurrency.pack(buffer, buffer_size, position, steps_per_sweep_);
  concurrency.pack(buffer, buffer_size, position, shifts_per_sweep_);
}

template <typename Concurrency>
void McSolverParameters<solver::SS_CT_HYB>::unpack(const Concurrency& concurrency, char* buffer,
                                                   int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, self_energy_tail_cutoff_);
  concurrency.unpack(buffer, buffer_size, position, steps_per_sweep_);
  concurrency.unpack(buffer, buffer_size, position, shifts_per_sweep_);
}

template <typename ReaderOrWriter>
void McSolverParameters<solver::SS_CT_HYB>::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("SS-CT-HYB");

    try {
      reader_or_writer.execute("self-energy-tail-cutoff", self_energy_tail_cutoff_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("steps-per-sweep", steps_per_sweep_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("shifts-per-sweep", shifts_per_sweep_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}

// Specialization for CT-INT
template <>
class McSolverParameters<solver::CT_INT> {
public:
  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  int getInitialConfigurationSize() const {
    return initial_configuration_size_;
  }

  std::array<double, 3> getAlphas() const {
    return std::array<double, 3>{alpha_dd_pos_, alpha_dd_neg_, alpha_ndd_};
  }

  double getDoubleUpdateProb() const {
    return double_update_prob_;
  }

  void setDoubleUpdateProb(const double p){
    double_update_prob_ = p;
  }

  bool doubleCountedInteraction() const {
    return double_counted_interaction_;
  }
  bool computeError() const {
    return compute_error_;
  }

private:
  int initial_configuration_size_ = 0;
  double alpha_dd_pos_ = 0.501;
  double alpha_dd_neg_ = 0;
  double alpha_ndd_ = 1e-4;
  double double_update_prob_= 0;
  bool double_counted_interaction_ = true;
  bool compute_error_ = false;
};

template <typename Concurrency>
int McSolverParameters<solver::CT_INT>::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;
  buffer_size += concurrency.get_buffer_size(initial_configuration_size_);
  buffer_size += concurrency.get_buffer_size(alpha_dd_pos_);
  buffer_size += concurrency.get_buffer_size(alpha_dd_neg_);
  buffer_size += concurrency.get_buffer_size(alpha_ndd_);
  buffer_size += concurrency.get_buffer_size(double_update_prob_);
  buffer_size += concurrency.get_buffer_size(double_counted_interaction_);
  buffer_size += concurrency.get_buffer_size(compute_error_);

  return buffer_size;
}

template <typename Concurrency>
void McSolverParameters<solver::CT_INT>::pack(const Concurrency& concurrency, char* buffer,
                                              const int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, initial_configuration_size_);
  concurrency.pack(buffer, buffer_size, position, alpha_dd_pos_);
  concurrency.pack(buffer, buffer_size, position, alpha_dd_neg_);
  concurrency.pack(buffer, buffer_size, position, alpha_ndd_);
  concurrency.pack(buffer, buffer_size, position, double_update_prob_);
  concurrency.pack(buffer, buffer_size, position, double_counted_interaction_);
  concurrency.pack(buffer, buffer_size, position, compute_error_);
}

template <typename Concurrency>
void McSolverParameters<solver::CT_INT>::unpack(const Concurrency& concurrency, char* buffer,
                                                const int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, initial_configuration_size_);
  concurrency.unpack(buffer, buffer_size, position, alpha_dd_pos_);
  concurrency.unpack(buffer, buffer_size, position, alpha_dd_neg_);
  concurrency.unpack(buffer, buffer_size, position, alpha_ndd_);
  concurrency.unpack(buffer, buffer_size, position, double_update_prob_);
  concurrency.unpack(buffer, buffer_size, position, double_counted_interaction_);
  concurrency.unpack(buffer, buffer_size, position, compute_error_);
}

// TODO: None of the open_group methods seem to throw.
template <typename ReaderOrWriter>
void McSolverParameters<solver::CT_INT>::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("CT-INT");

    auto tryToRead = [&](const std::string& name, auto& obj) {
      try {
        reader_or_writer.execute(name, obj);
      }
      catch (const std::exception& /*err*/) {
      }
    };

    tryToRead("initial-configuration-size", initial_configuration_size_);
    tryToRead("alpha-dd-pos", alpha_dd_pos_);
    tryToRead("alpha-dd-neg", alpha_dd_neg_);
    tryToRead("alpha-ndd", alpha_ndd_);
    tryToRead("double-update-probability", double_update_prob_);
    tryToRead("double-counted-interaction", double_counted_interaction_);
    tryToRead("compute-error-bar", compute_error_);
    reader_or_writer.close_group();
  }
  catch (const std::exception& /*r_e*/) {
    std::cout << "\nNo CT-INT solver parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  if(double_update_prob_ < 0 or double_update_prob_> 1.)
    throw(std::out_of_range("double-update-probability must be in [0,1]."));
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_MC_SOLVER_PARAMETERS_HPP
