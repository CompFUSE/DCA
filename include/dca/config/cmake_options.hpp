// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class stores and prints a list of DCA++'s CMake options and their values.

#ifndef DCA_CONFIG_CMAKE_OPTIONS_HPP
#define DCA_CONFIG_CMAKE_OPTIONS_HPP

#include <string>

namespace dca {
namespace config {
// dca::config::

struct CMakeOptions {
  // CUDA
  static const std::string dca_with_cuda;
  static const std::string cuda_gpu_arch;

  // Parallelization
  static const std::string dca_with_mpi;
  static const std::string dca_with_threaded_solver;
  static const std::string dca_threading_library;

  // Others
  static const std::string dca_cluster_solver;
  static const std::string dca_lattice;
  static const std::string dca_point_group;
  static const std::string dca_model;
  static const std::string dca_rng;
  // static const std::string dca_rng_class;
  // static const std::string dca_rng_header;
  // static const std::string dca_rng_library;
  static const std::string dca_profiler;
  static const std::string dca_with_autotuning;
  static const std::string dca_with_gnuplot;
  static const std::string dca_with_single_precision_mc;
  static const std::string dca_with_single_precision_tp_measurements;
  static const std::string dca_with_memory_savings;
  static const std::string dca_with_managed_memory;

  static const std::string dca_with_qmc_bit;

  static void print();
};

}  // namespace config
}  // namespace dca

#endif  // DCA_CONFIG_CMAKE_OPTIONS_HPP
