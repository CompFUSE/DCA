// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the DCA parameters.

#ifndef DCA_PHYS_PARAMETERS_DCA_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_DCA_PARAMETERS_HPP

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class DcaParameters {
public:
  DcaParameters()
      : initial_self_energy_("zero"),
        dca_iterations_(1),
        dca_accuracy_(0.),
        self_energy_mixing_factor_(1.),
        interacting_orbitals_{0},

        do_finite_size_qmc_(false),

        do_simple_q_points_summation_(true),
        k_mesh_recursion_(0),
        coarsegraining_periods_(0),
        quadrature_rule_(1),
        coarsegraining_threads_(1),
        tail_frequencies_(0),

        do_dca_plus_(false),
        deconvolution_iterations_(16),
        deconvolution_tolerance_(1.e-3),
        hts_approximation_(false),
        hts_threads_(1) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const std::string& get_initial_self_energy() const {
    return initial_self_energy_;
  }
  int get_dca_iterations() const {
    return dca_iterations_;
  }
  void set_dca_iterations(const int iterations) {
    assert(iterations > 0);
    dca_iterations_ = iterations;
  }
  double get_dca_accuracy() const {
    return dca_accuracy_;
  }
  double get_self_energy_mixing_factor() const {
    return self_energy_mixing_factor_;
  }
  const std::vector<int>& get_interacting_orbitals() const {
    return interacting_orbitals_;
  }
  bool do_finite_size_qmc() const {
    return do_finite_size_qmc_;
  }

  // Compute the coarse grained G function as a simple average over the q points. This will be
  // computationally faster and perform better on multiband systems.
  // Currently not supported for DCA+.
  bool do_simple_q_points_summation() const {
    return do_simple_q_points_summation_;
  }
  int get_k_mesh_recursion() const {
    return k_mesh_recursion_;
  }
  int get_coarsegraining_periods() const {
    return coarsegraining_periods_;
  }
  int get_quadrature_rule() const {
    return quadrature_rule_;
  }
  int get_coarsegraining_threads() const {
    return coarsegraining_threads_;
  }
  int get_tail_frequencies() const {
    return tail_frequencies_;
  }
  bool do_dca_plus() const {
    return do_dca_plus_;
  }
  int get_deconvolution_iterations() const {
    return deconvolution_iterations_;
  }
  double get_deconvolution_tolerance() const {
    return deconvolution_tolerance_;
  }
  bool hts_approximation() const {
    return hts_approximation_;
  }
  int get_hts_threads() const {
    return hts_threads_;
  }

private:
  std::string initial_self_energy_;
  int dca_iterations_;
  double dca_accuracy_;
  double self_energy_mixing_factor_;
  std::vector<int> interacting_orbitals_;

  bool do_finite_size_qmc_;

  // coarse-graining
  bool do_simple_q_points_summation_;
  int k_mesh_recursion_;
  int coarsegraining_periods_;
  int quadrature_rule_;
  int coarsegraining_threads_;
  int tail_frequencies_;

  // DCA+
  bool do_dca_plus_;
  int deconvolution_iterations_;
  double deconvolution_tolerance_;
  bool hts_approximation_;
  int hts_threads_;
};

template <typename Concurrency>
int DcaParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(initial_self_energy_);
  buffer_size += concurrency.get_buffer_size(dca_iterations_);
  buffer_size += concurrency.get_buffer_size(dca_accuracy_);
  buffer_size += concurrency.get_buffer_size(self_energy_mixing_factor_);
  buffer_size += concurrency.get_buffer_size(interacting_orbitals_);
  buffer_size += concurrency.get_buffer_size(do_finite_size_qmc_);
  buffer_size += concurrency.get_buffer_size(do_simple_q_points_summation_);
  buffer_size += concurrency.get_buffer_size(k_mesh_recursion_);
  buffer_size += concurrency.get_buffer_size(coarsegraining_periods_);
  buffer_size += concurrency.get_buffer_size(quadrature_rule_);
  buffer_size += concurrency.get_buffer_size(coarsegraining_threads_);
  buffer_size += concurrency.get_buffer_size(tail_frequencies_);
  buffer_size += concurrency.get_buffer_size(do_dca_plus_);
  buffer_size += concurrency.get_buffer_size(deconvolution_iterations_);
  buffer_size += concurrency.get_buffer_size(deconvolution_tolerance_);
  buffer_size += concurrency.get_buffer_size(hts_approximation_);
  buffer_size += concurrency.get_buffer_size(hts_threads_);

  return buffer_size;
}

template <typename Concurrency>
void DcaParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                         int& position) const {
  concurrency.pack(buffer, buffer_size, position, initial_self_energy_);
  concurrency.pack(buffer, buffer_size, position, dca_iterations_);
  concurrency.pack(buffer, buffer_size, position, dca_accuracy_);
  concurrency.pack(buffer, buffer_size, position, self_energy_mixing_factor_);
  concurrency.pack(buffer, buffer_size, position, interacting_orbitals_);
  concurrency.pack(buffer, buffer_size, position, do_finite_size_qmc_);
  concurrency.pack(buffer, buffer_size, position, do_simple_q_points_summation_);
  concurrency.pack(buffer, buffer_size, position, k_mesh_recursion_);
  concurrency.pack(buffer, buffer_size, position, coarsegraining_periods_);
  concurrency.pack(buffer, buffer_size, position, quadrature_rule_);
  concurrency.pack(buffer, buffer_size, position, coarsegraining_threads_);
  concurrency.pack(buffer, buffer_size, position, tail_frequencies_);
  concurrency.pack(buffer, buffer_size, position, do_dca_plus_);
  concurrency.pack(buffer, buffer_size, position, deconvolution_iterations_);
  concurrency.pack(buffer, buffer_size, position, deconvolution_tolerance_);
  concurrency.pack(buffer, buffer_size, position, hts_approximation_);
  concurrency.pack(buffer, buffer_size, position, hts_threads_);
}

template <typename Concurrency>
void DcaParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                           int& position) {
  concurrency.unpack(buffer, buffer_size, position, initial_self_energy_);
  concurrency.unpack(buffer, buffer_size, position, dca_iterations_);
  concurrency.unpack(buffer, buffer_size, position, dca_accuracy_);
  concurrency.unpack(buffer, buffer_size, position, self_energy_mixing_factor_);
  concurrency.unpack(buffer, buffer_size, position, interacting_orbitals_);
  concurrency.unpack(buffer, buffer_size, position, do_finite_size_qmc_);
  concurrency.unpack(buffer, buffer_size, position, do_simple_q_points_summation_);
  concurrency.unpack(buffer, buffer_size, position, k_mesh_recursion_);
  concurrency.unpack(buffer, buffer_size, position, coarsegraining_periods_);
  concurrency.unpack(buffer, buffer_size, position, quadrature_rule_);
  concurrency.unpack(buffer, buffer_size, position, coarsegraining_threads_);
  concurrency.unpack(buffer, buffer_size, position, tail_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, do_dca_plus_);
  concurrency.unpack(buffer, buffer_size, position, deconvolution_iterations_);
  concurrency.unpack(buffer, buffer_size, position, deconvolution_tolerance_);
  concurrency.unpack(buffer, buffer_size, position, hts_approximation_);
  concurrency.unpack(buffer, buffer_size, position, hts_threads_);
}

template <typename ReaderOrWriter>
void DcaParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  auto try_to_read = [&](const std::string& name, auto& var) {
    try {
      reader_or_writer.execute(name, var);
    }
    catch (const std::exception& r_e) {
    }
  };

  try {
    reader_or_writer.open_group("DCA");

    try_to_read("initial-self-energy", initial_self_energy_);
    try_to_read("iterations", dca_iterations_);
    try_to_read("accuracy", dca_accuracy_);
    try_to_read("self-energy-mixing-factor", self_energy_mixing_factor_);
    try_to_read("interacting-orbitals", interacting_orbitals_);

    try_to_read("do-finite-size-QMC", do_finite_size_qmc_);

    try {
      reader_or_writer.open_group("coarse-graining");

      try_to_read("do-simple-q-points-summation", do_simple_q_points_summation_);
      try_to_read("k-mesh-recursion", k_mesh_recursion_);
      try_to_read("periods", coarsegraining_periods_);
      try_to_read("quadrature-rule", quadrature_rule_);
      try_to_read("threads", coarsegraining_threads_);
      try_to_read("tail-frequencies", tail_frequencies_);

      reader_or_writer.close_group();
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.open_group("DCA+");

      try_to_read("do-DCA+", do_dca_plus_);
      try_to_read("deconvolution-iterations", deconvolution_iterations_);
      try_to_read("deconvolution-tolerance", deconvolution_tolerance_);
      try_to_read("HTS-approximation", hts_approximation_);
      try_to_read("HTS-threads", hts_threads_);

      reader_or_writer.close_group();
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }

  // Check read parameters for consistency.
  if (reader_or_writer.is_reader()) {
    if (do_finite_size_qmc_ && do_dca_plus_)
      throw std::logic_error("Finite-size QMC and DCA+ are mutually exclusive options.");
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_DCA_PARAMETERS_HPP
