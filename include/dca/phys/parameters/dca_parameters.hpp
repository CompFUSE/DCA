// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the DCA parameters.

#ifndef DCA_PHYS_PARAMETERS_DCA_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_DCA_PARAMETERS_HPP

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
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const std::string& get_initial_self_energy() const {
    return initial_self_energy_;
  }
  int get_dca_iterations() const {
    return dca_iterations_;
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
void DcaParameters::pack(const Concurrency& concurrency, int* buffer, int buffer_size,
                         int& position) const {
  concurrency.pack(buffer, buffer_size, position, initial_self_energy_);
  concurrency.pack(buffer, buffer_size, position, dca_iterations_);
  concurrency.pack(buffer, buffer_size, position, dca_accuracy_);
  concurrency.pack(buffer, buffer_size, position, self_energy_mixing_factor_);
  concurrency.pack(buffer, buffer_size, position, interacting_orbitals_);
  concurrency.pack(buffer, buffer_size, position, do_finite_size_qmc_);
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
void DcaParameters::unpack(const Concurrency& concurrency, int* buffer, int buffer_size,
                           int& position) {
  concurrency.unpack(buffer, buffer_size, position, initial_self_energy_);
  concurrency.unpack(buffer, buffer_size, position, dca_iterations_);
  concurrency.unpack(buffer, buffer_size, position, dca_accuracy_);
  concurrency.unpack(buffer, buffer_size, position, self_energy_mixing_factor_);
  concurrency.unpack(buffer, buffer_size, position, interacting_orbitals_);
  concurrency.unpack(buffer, buffer_size, position, do_finite_size_qmc_);
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
  try {
    reader_or_writer.open_group("DCA");

    try {
      reader_or_writer.execute("initial-self-energy", initial_self_energy_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("iterations", dca_iterations_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("accuracy", dca_accuracy_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("self-energy-mixing-factor", self_energy_mixing_factor_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("interacting-orbitals", interacting_orbitals_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("do-finite-size-QMC", do_finite_size_qmc_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.open_group("coarse-graining");

      try {
        reader_or_writer.execute("k-mesh-recursion", k_mesh_recursion_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("periods", coarsegraining_periods_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("quadrature-rule", quadrature_rule_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("threads", coarsegraining_threads_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("tail-frequencies", tail_frequencies_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.open_group("DCA+");

      try {
        reader_or_writer.execute("do-DCA+", do_dca_plus_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("deconvolution-iterations", deconvolution_iterations_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("deconvolution-tolerance", deconvolution_tolerance_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("HTS-approximation", hts_approximation_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("HTS-threads", hts_threads_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
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

#endif  // DCA_PHYS_PARAMETERS_DCA_PARAMETERS_HPP
