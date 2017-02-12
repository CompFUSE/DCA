// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class reads, stores, and writes the Monte Carlo solver parameters.
// It is templated on the MC solver name and only implemented when specialized for CT-AUX and
// SS-CT-HYB.

#ifndef DCA_PHYS_PARAMETERS_MC_SOLVER_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_MC_SOLVER_PARAMETERS_HPP

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
        initial_matrix_size_(128),
        submatrix_size_(128),
        additional_time_measurements_(false) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_expansion_parameter_K() const {
    return expansion_parameter_K_;
  }
  int get_initial_matrix_size() const {
    return initial_matrix_size_;
  }
  int get_submatrix_size() const {
    return submatrix_size_;
  }
  bool additional_time_measurements() const {
    return additional_time_measurements_;
  }

private:
  double expansion_parameter_K_;
  int initial_matrix_size_;
  int submatrix_size_;
  bool additional_time_measurements_;
};

template <typename Concurrency>
int McSolverParameters<solver::CT_AUX>::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(expansion_parameter_K_);
  buffer_size += concurrency.get_buffer_size(initial_matrix_size_);
  buffer_size += concurrency.get_buffer_size(submatrix_size_);
  buffer_size += concurrency.get_buffer_size(additional_time_measurements_);

  return buffer_size;
}

template <typename Concurrency>
void McSolverParameters<solver::CT_AUX>::pack(const Concurrency& concurrency, int* buffer,
                                              int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, expansion_parameter_K_);
  concurrency.pack(buffer, buffer_size, position, initial_matrix_size_);
  concurrency.pack(buffer, buffer_size, position, submatrix_size_);
  concurrency.pack(buffer, buffer_size, position, additional_time_measurements_);
}

template <typename Concurrency>
void McSolverParameters<solver::CT_AUX>::unpack(const Concurrency& concurrency, int* buffer,
                                                int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, expansion_parameter_K_);
  concurrency.unpack(buffer, buffer_size, position, initial_matrix_size_);
  concurrency.unpack(buffer, buffer_size, position, submatrix_size_);
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
      reader_or_writer.execute("initial-matrix-size", initial_matrix_size_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("submatrix-size", submatrix_size_);
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
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

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
void McSolverParameters<solver::SS_CT_HYB>::pack(const Concurrency& concurrency, int* buffer,
                                                 int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, self_energy_tail_cutoff_);
  concurrency.pack(buffer, buffer_size, position, steps_per_sweep_);
  concurrency.pack(buffer, buffer_size, position, shifts_per_sweep_);
}

template <typename Concurrency>
void McSolverParameters<solver::SS_CT_HYB>::unpack(const Concurrency& concurrency, int* buffer,
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

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_MC_SOLVER_PARAMETERS_HPP
