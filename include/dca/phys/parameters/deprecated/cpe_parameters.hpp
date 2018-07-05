// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the Continuous-pole-expansion (CPE) parameters.

#ifndef DCA_PHYS_PARAMETERS_CPE_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_CPE_PARAMETERS_HPP

#include <stdexcept>
#include <string>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class CpeParameters {
public:
  CpeParameters()
      : do_CPE_("false"),
        N_wn_(64),
        CPE_smoothing_factor_(1.),
        max_CPE_iterations_(100),
        max_CPE_error_(0.),
        simulate_gaussian_noise_("false"),
        nr_of_CPE_samples_(1),
        simulated_CPE_stddev_(0.),
        compute_free_spectrum_("false"),
        compute_lattice_spectrum_("false"),
        compute_cluster_spectrum_("false") {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  bool do_CPE() const {
    return (do_CPE_ == "true");
  }
  int get_N_wn() {
    return N_wn_;
  }
  double get_CPE_smoothing_factor() const {
    return CPE_smoothing_factor_;
  }
  int get_max_CPE_iterations() const {
    return max_CPE_iterations_;
  }
  double get_max_CPE_error() const {
    return max_CPE_error_;
  }
  bool simulate_gaussian_noise() const {
    return (simulate_gaussian_noise_ == "true");
  }
  int get_nr_of_CPE_samples() const {
    return nr_of_CPE_samples_;
  }
  double get_simulated_CPE_stddev() const {
    return simulated_CPE_stddev_;
  }
  bool compute_free_spectrum() const {
    return (compute_free_spectrum_ == "true");
  }
  bool compute_lattice_spectrum() const {
    return (compute_lattice_spectrum_ == "true");
  }
  bool compute_cluster_spectrum() const {
    return (compute_cluster_spectrum_ == "true");
  }

private:
  std::string do_CPE_;
  int N_wn_;
  double CPE_smoothing_factor_;
  int max_CPE_iterations_;
  double max_CPE_error_;
  std::string simulate_gaussian_noise_;
  int nr_of_CPE_samples_;
  double simulated_CPE_stddev_;
  std::string compute_free_spectrum_;
  std::string compute_lattice_spectrum_;
  std::string compute_cluster_spectrum_;
};

template <typename Concurrency>
int CpeParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(do_CPE_);
  buffer_size += concurrency.get_buffer_size(N_wn_);
  buffer_size += concurrency.get_buffer_size(CPE_smoothing_factor_);
  buffer_size += concurrency.get_buffer_size(max_CPE_iterations_);
  buffer_size += concurrency.get_buffer_size(max_CPE_error_);
  buffer_size += concurrency.get_buffer_size(simulate_gaussian_noise_);
  buffer_size += concurrency.get_buffer_size(nr_of_CPE_samples_);
  buffer_size += concurrency.get_buffer_size(simulated_CPE_stddev_);
  buffer_size += concurrency.get_buffer_size(compute_free_spectrum_);
  buffer_size += concurrency.get_buffer_size(compute_lattice_spectrum_);
  buffer_size += concurrency.get_buffer_size(compute_cluster_spectrum_);

  return buffer_size;
}

template <typename Concurrency>
void CpeParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                         int& position) const {
  concurrency.pack(buffer, buffer_size, position, do_CPE_);
  concurrency.pack(buffer, buffer_size, position, N_wn_);
  concurrency.pack(buffer, buffer_size, position, CPE_smoothing_factor_);
  concurrency.pack(buffer, buffer_size, position, max_CPE_iterations_);
  concurrency.pack(buffer, buffer_size, position, max_CPE_error_);
  concurrency.pack(buffer, buffer_size, position, simulate_gaussian_noise_);
  concurrency.pack(buffer, buffer_size, position, nr_of_CPE_samples_);
  concurrency.pack(buffer, buffer_size, position, simulated_CPE_stddev_);
  concurrency.pack(buffer, buffer_size, position, compute_free_spectrum_);
  concurrency.pack(buffer, buffer_size, position, compute_lattice_spectrum_);
  concurrency.pack(buffer, buffer_size, position, compute_cluster_spectrum_);
}

template <typename Concurrency>
void CpeParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                           int& position) {
  concurrency.unpack(buffer, buffer_size, position, do_CPE_);
  concurrency.unpack(buffer, buffer_size, position, N_wn_);
  concurrency.unpack(buffer, buffer_size, position, CPE_smoothing_factor_);
  concurrency.unpack(buffer, buffer_size, position, max_CPE_iterations_);
  concurrency.unpack(buffer, buffer_size, position, max_CPE_error_);
  concurrency.unpack(buffer, buffer_size, position, simulate_gaussian_noise_);
  concurrency.unpack(buffer, buffer_size, position, nr_of_CPE_samples_);
  concurrency.unpack(buffer, buffer_size, position, simulated_CPE_stddev_);
  concurrency.unpack(buffer, buffer_size, position, compute_free_spectrum_);
  concurrency.unpack(buffer, buffer_size, position, compute_lattice_spectrum_);
  concurrency.unpack(buffer, buffer_size, position, compute_cluster_spectrum_);
}

template <typename ReaderOrWriter>
void CpeParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("CPE-parameters");

    try {
      reader_or_writer.execute("do-CPE", do_CPE_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("number-of-matsubara-freqencies", N_wn_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("smoothing-factor", CPE_smoothing_factor_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("max-CPE-iterations", max_CPE_iterations_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("max-CPE-error", max_CPE_error_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("simulate-Gaussian-noise", simulate_gaussian_noise_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("nr-of-samples", nr_of_CPE_samples_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("simulated-stddev", simulated_CPE_stddev_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("compute-free-spectrum", compute_free_spectrum_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("compute-lattice-spectrum", compute_lattice_spectrum_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("compute-cluster-spectrum", compute_cluster_spectrum_);
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

#endif  // DCA_PHYS_PARAMETERS_CPE_PARAMETERS_HPP
