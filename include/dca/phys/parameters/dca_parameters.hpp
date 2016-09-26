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
  DcaParameters(int lattice_dimension)
      : do_DCA_plus_("false"),
        interacting_bands_(0),
        DCA_iterations_(1),
        DCA_accuracy_(1.e-5),
        DCA_mixing_factor_(1.),
        DCA_cluster_(lattice_dimension, std::vector<int>(lattice_dimension, 0)),

        k_mesh_refinement_(3),
        number_of_periods_(0),
        quadrature_rule_(3),
        precompute_Hamiltonian_("true"),
        nr_coarsegraining_threads_(1),
        nr_tail_frequencies_(0),
        phi_k_integration_accuracy_(1.e-3),
        print_phi_k_("false"),

        interpolation_method_("wannier-interpolation"),
        HTS_approximation_("false"),
        deconvolution_tolerance_(1.e-2),
        max_deconvolution_iterations_(16) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  bool do_DCA_plus() const {
    return (do_DCA_plus_ == "true");
  }
  const std::vector<int>& get_interacting_bands() const {
    return interacting_bands_;
  }
  bool is_an_interacting_band(int i) const {
    for (std::size_t l = 0; l < interacting_bands_.size(); l++) {
      if (i == interacting_bands_[l])
        return true;
    }
    return false;
  }
  int get_DCA_iterations() const {
    return DCA_iterations_;
  }
  double get_DCA_accuracy() const {
    return DCA_accuracy_;
  }
  double get_DCA_mixing_factor() const {
    return DCA_mixing_factor_;
  }
  const std::vector<std::vector<int>>& get_DCA_cluster() const {
    return DCA_cluster_;
  }

  int get_k_mesh_refinement() const {
    return k_mesh_refinement_;
  }
  int get_number_of_periods() const {
    return number_of_periods_;
  }
  int get_quadrature_rule() const {
    return quadrature_rule_;
  }
  bool precompute_Hamiltonian() const {
    return (precompute_Hamiltonian_ == "true");
  }
  int get_nr_coarsegraining_threads() const {
    return nr_coarsegraining_threads_;
  }
  int get_number_of_tail_frequencies() const {
    return nr_tail_frequencies_;
  }
  double get_phi_k_integration_accuracy() const {
    return phi_k_integration_accuracy_;
  }
  bool print_phi_k() const {
    return (print_phi_k_ == "true");
  }

  const std::string& get_interpolation_method() const {
    return interpolation_method_;
  }
  bool use_HTS_approximation() const {
    return (HTS_approximation_ == "true");
  }
  double get_deconvolution_tolerance() const {
    return deconvolution_tolerance_;
  }
  int get_max_deconvolution_iterations() const {
    return max_deconvolution_iterations_;
  }

private:
  // general
  std::string do_DCA_plus_;
  std::vector<int> interacting_bands_;
  int DCA_iterations_;
  double DCA_accuracy_;
  double DCA_mixing_factor_;
  std::vector<std::vector<int>> DCA_cluster_;

  // cluster-mapping
  int k_mesh_refinement_;
  int number_of_periods_;
  int quadrature_rule_;
  std::string precompute_Hamiltonian_;
  int nr_coarsegraining_threads_;
  int nr_tail_frequencies_;
  double phi_k_integration_accuracy_;
  std::string print_phi_k_;

  // lattice-mapping
  std::string interpolation_method_;
  std::string HTS_approximation_;
  double deconvolution_tolerance_;
  int max_deconvolution_iterations_;
};

template <typename Concurrency>
int DcaParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(do_DCA_plus_);
  buffer_size += concurrency.get_buffer_size(interacting_bands_);
  buffer_size += concurrency.get_buffer_size(DCA_iterations_);
  buffer_size += concurrency.get_buffer_size(DCA_accuracy_);
  buffer_size += concurrency.get_buffer_size(DCA_mixing_factor_);
  buffer_size += concurrency.get_buffer_size(DCA_cluster_);

  buffer_size += concurrency.get_buffer_size(k_mesh_refinement_);
  buffer_size += concurrency.get_buffer_size(number_of_periods_);
  buffer_size += concurrency.get_buffer_size(quadrature_rule_);
  buffer_size += concurrency.get_buffer_size(precompute_Hamiltonian_);
  buffer_size += concurrency.get_buffer_size(nr_coarsegraining_threads_);
  buffer_size += concurrency.get_buffer_size(nr_tail_frequencies_);
  buffer_size += concurrency.get_buffer_size(phi_k_integration_accuracy_);
  buffer_size += concurrency.get_buffer_size(print_phi_k_);

  buffer_size += concurrency.get_buffer_size(interpolation_method_);
  buffer_size += concurrency.get_buffer_size(HTS_approximation_);
  buffer_size += concurrency.get_buffer_size(deconvolution_tolerance_);
  buffer_size += concurrency.get_buffer_size(max_deconvolution_iterations_);

  return buffer_size;
}

template <typename Concurrency>
void DcaParameters::pack(const Concurrency& concurrency, int* buffer, int buffer_size,
                         int& position) const {
  concurrency.pack(buffer, buffer_size, position, do_DCA_plus_);
  concurrency.pack(buffer, buffer_size, position, interacting_bands_);
  concurrency.pack(buffer, buffer_size, position, DCA_iterations_);
  concurrency.pack(buffer, buffer_size, position, DCA_accuracy_);
  concurrency.pack(buffer, buffer_size, position, DCA_mixing_factor_);
  concurrency.pack(buffer, buffer_size, position, DCA_cluster_);

  concurrency.pack(buffer, buffer_size, position, k_mesh_refinement_);
  concurrency.pack(buffer, buffer_size, position, number_of_periods_);
  concurrency.pack(buffer, buffer_size, position, quadrature_rule_);
  concurrency.pack(buffer, buffer_size, position, precompute_Hamiltonian_);
  concurrency.pack(buffer, buffer_size, position, nr_coarsegraining_threads_);
  concurrency.pack(buffer, buffer_size, position, nr_tail_frequencies_);
  concurrency.pack(buffer, buffer_size, position, phi_k_integration_accuracy_);
  concurrency.pack(buffer, buffer_size, position, print_phi_k_);

  concurrency.pack(buffer, buffer_size, position, interpolation_method_);
  concurrency.pack(buffer, buffer_size, position, HTS_approximation_);
  concurrency.pack(buffer, buffer_size, position, deconvolution_tolerance_);
  concurrency.pack(buffer, buffer_size, position, max_deconvolution_iterations_);
}

template <typename Concurrency>
void DcaParameters::unpack(const Concurrency& concurrency, int* buffer, int buffer_size,
                           int& position) {
  concurrency.unpack(buffer, buffer_size, position, do_DCA_plus_);
  concurrency.unpack(buffer, buffer_size, position, interacting_bands_);
  concurrency.unpack(buffer, buffer_size, position, DCA_iterations_);
  concurrency.unpack(buffer, buffer_size, position, DCA_accuracy_);
  concurrency.unpack(buffer, buffer_size, position, DCA_mixing_factor_);
  concurrency.unpack(buffer, buffer_size, position, DCA_cluster_);

  concurrency.unpack(buffer, buffer_size, position, k_mesh_refinement_);
  concurrency.unpack(buffer, buffer_size, position, number_of_periods_);
  concurrency.unpack(buffer, buffer_size, position, quadrature_rule_);
  concurrency.unpack(buffer, buffer_size, position, precompute_Hamiltonian_);
  concurrency.unpack(buffer, buffer_size, position, nr_coarsegraining_threads_);
  concurrency.unpack(buffer, buffer_size, position, nr_tail_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, phi_k_integration_accuracy_);
  concurrency.unpack(buffer, buffer_size, position, print_phi_k_);

  concurrency.unpack(buffer, buffer_size, position, interpolation_method_);
  concurrency.unpack(buffer, buffer_size, position, HTS_approximation_);
  concurrency.unpack(buffer, buffer_size, position, deconvolution_tolerance_);
  concurrency.unpack(buffer, buffer_size, position, max_deconvolution_iterations_);
}

template <typename ReaderOrWriter>
void DcaParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("DCA");

    try {
      reader_or_writer.execute("do-DCA+", do_DCA_plus_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("interacting-bands", interacting_bands_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("DCA-iterations", DCA_iterations_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("DCA-accuracy", DCA_accuracy_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("DCA-mixing-factor", DCA_mixing_factor_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("cluster", DCA_cluster_);
    }
    catch (const std::exception& r_e) {
    }

    {
      reader_or_writer.open_group("cluster-mapping");

      try {
        reader_or_writer.execute("k-mesh-refinement", k_mesh_refinement_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("number-of-periods", number_of_periods_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("quadrature-rule", quadrature_rule_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("precompute-Hamiltonian", precompute_Hamiltonian_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("number-of-threads", nr_coarsegraining_threads_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("number-of-tail-frequencies", nr_tail_frequencies_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("phi(k) integration accuracy", phi_k_integration_accuracy_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("print-phi(k)", print_phi_k_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }

    {
      reader_or_writer.open_group("lattice-mapping");

      try {
        reader_or_writer.execute("interpolation-method", interpolation_method_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("HTS-approximation", HTS_approximation_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("deconvolution-tolerance", deconvolution_tolerance_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("max-deconvolution-iterations", max_deconvolution_iterations_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\nNo DCA parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_DCA_PARAMETERS_HPP
