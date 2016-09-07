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
//
// TODO: Const correctness.

#ifndef DCA_PHYS_PARAMETERS_MC_SOLVER_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_MC_SOLVER_PARAMETERS_HPP

#include <iostream>
#include <stdexcept>
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_name.hpp"

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

// Empty class template
template <DCA::ClusterSolverName solver_name>
class McSolverParameters {};

// Specialization for CT-AUX
template <>
class McSolverParameters<DCA::CT_AUX_CLUSTER_SOLVER> {
public:
  McSolverParameters() : submatrix_size_(128), initial_matrix_size_(128), K_parameter_(1.) {}

  template <typename Concurrency>
  int getBufferSize(Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  int get_submatrix_size() const {
    return submatrix_size_;
  }
  int get_initial_matrix_size() const {
    return initial_matrix_size_;
  }
  double get_K_parameter() const {
    return K_parameter_;
  }

private:
  int submatrix_size_;
  int initial_matrix_size_;
  double K_parameter_;
};

template <typename Concurrency>
int McSolverParameters<DCA::CT_AUX_CLUSTER_SOLVER>::getBufferSize(Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(submatrix_size_);
  buffer_size += concurrency.get_buffer_size(initial_matrix_size_);
  buffer_size += concurrency.get_buffer_size(K_parameter_);

  return buffer_size;
}

template <typename Concurrency>
void McSolverParameters<DCA::CT_AUX_CLUSTER_SOLVER>::pack(Concurrency& concurrency, int* buffer,
                                                          int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, submatrix_size_);
  concurrency.pack(buffer, buffer_size, position, initial_matrix_size_);
  concurrency.pack(buffer, buffer_size, position, K_parameter_);
}

template <typename Concurrency>
void McSolverParameters<DCA::CT_AUX_CLUSTER_SOLVER>::unpack(Concurrency& concurrency, int* buffer,
                                                            int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, submatrix_size_);
  concurrency.unpack(buffer, buffer_size, position, initial_matrix_size_);
  concurrency.unpack(buffer, buffer_size, position, K_parameter_);
}

// TODO: None of the open_group methods seem to throw.
template <typename ReaderOrWriter>
void McSolverParameters<DCA::CT_AUX_CLUSTER_SOLVER>::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("CT-AUX-solver");

    try {
      reader_or_writer.execute("submatrix-size", submatrix_size_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("initial-matrix-size", initial_matrix_size_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("K-parameter", K_parameter_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\nNo CT-AUX solver parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

// Specialization for SS-CT-HYB
template <>
class McSolverParameters<DCA::SS_CT_HYB> {
public:
  McSolverParameters()
      : Sigma_tail_cutoff_(0), steps_per_sweep_(0.4), swaps_per_sweep_(0.1), shifts_per_sweep_(0.4) {}

  template <typename Concurrency>
  int getBufferSize(Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  int get_Sigma_tail_cutoff() const {
    return Sigma_tail_cutoff_;
  }
  double get_steps_per_sweep() const {
    return steps_per_sweep_;
  }
  double get_swaps_per_sweep() const {
    return swaps_per_sweep_;
  }
  double get_shifts_per_sweep() const {
    return shifts_per_sweep_;
  }

private:
  int Sigma_tail_cutoff_;
  double steps_per_sweep_;
  double swaps_per_sweep_;
  double shifts_per_sweep_;
};

template <typename Concurrency>
int McSolverParameters<DCA::SS_CT_HYB>::getBufferSize(Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(Sigma_tail_cutoff_);
  buffer_size += concurrency.get_buffer_size(steps_per_sweep_);
  buffer_size += concurrency.get_buffer_size(swaps_per_sweep_);
  buffer_size += concurrency.get_buffer_size(shifts_per_sweep_);

  return buffer_size;
}

template <typename Concurrency>
void McSolverParameters<DCA::SS_CT_HYB>::pack(Concurrency& concurrency, int* buffer,
                                              int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, Sigma_tail_cutoff_);
  concurrency.pack(buffer, buffer_size, position, steps_per_sweep_);
  concurrency.pack(buffer, buffer_size, position, swaps_per_sweep_);
  concurrency.pack(buffer, buffer_size, position, shifts_per_sweep_);
}

template <typename Concurrency>
void McSolverParameters<DCA::SS_CT_HYB>::unpack(Concurrency& concurrency, int* buffer,
                                                int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, Sigma_tail_cutoff_);
  concurrency.unpack(buffer, buffer_size, position, steps_per_sweep_);
  concurrency.unpack(buffer, buffer_size, position, swaps_per_sweep_);
  concurrency.unpack(buffer, buffer_size, position, shifts_per_sweep_);
}

template <typename ReaderOrWriter>
void McSolverParameters<DCA::SS_CT_HYB>::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("SS-CT-HYB-solver");

    try {
      reader_or_writer.execute("Sigma-tail-cutoff", Sigma_tail_cutoff_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("steps-per-sweep", steps_per_sweep_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("swaps-per-sweep", swaps_per_sweep_);
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
    std::cout << "\nNo SS-CT-HYB solver parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_MC_SOLVER_PARAMETERS_HPP
