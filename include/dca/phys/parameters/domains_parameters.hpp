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
// This class reads, stores, and writes parameters that determine the domains.

#ifndef DCA_PHYS_PARAMETERS_DOMAINS_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_DOMAINS_PARAMETERS_HPP

#include <vector>
#include <stdexcept>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class DomainsParameters {
public:
  DomainsParameters(int dimension)
      : cluster_(dimension, std::vector<int>(dimension, 0)),
        sp_host_(dimension, std::vector<int>(dimension, 0)),
        tp_host_(dimension, std::vector<int>(dimension, 0)),
        sp_time_intervals_(128),
        time_intervals_for_time_measurements_(1),
        sp_fermionic_frequencies_(256),
        hts_bosonic_frequencies_(0),
        four_point_fermionic_frequencies_(1),
        min_real_frequency_(-10.),
        max_real_frequency_(10.),
        real_frequencies_(3),
        imaginary_damping_(0.01) {
    // Set all grids to their default value, the lattice basis.
    for (int i = 0; i < dimension; ++i) {
      cluster_[i][i] = 1;
      sp_host_[i][i] = 1;
      tp_host_[i][i] = 1;
    }
  }

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const std::vector<std::vector<int>>& get_cluster() const {
    return cluster_;
  }
  const std::vector<std::vector<int>>& get_sp_host() const {
    return sp_host_;
  }
  const std::vector<std::vector<int>>& get_tp_host() const {
    return tp_host_;
  }
  int get_sp_time_intervals() const {
    return sp_time_intervals_;
  }
  int get_time_intervals_for_time_measurements() const {
    return time_intervals_for_time_measurements_;
  }
  int get_sp_fermionic_frequencies() const {
    return sp_fermionic_frequencies_;
  }
  int get_hts_bosonic_frequencies() const {
    return hts_bosonic_frequencies_;
  }
  int get_four_point_fermionic_frequencies() const {
    return four_point_fermionic_frequencies_;
  }
  double get_min_real_frequency() const {
    return min_real_frequency_;
  }
  double get_max_real_frequency() const {
    return max_real_frequency_;
  }
  int get_real_frequencies() const {
    return real_frequencies_;
  }
  double get_imaginary_damping() const {
    return imaginary_damping_;
  }

private:
  std::vector<std::vector<int>> cluster_;
  std::vector<std::vector<int>> sp_host_;
  std::vector<std::vector<int>> tp_host_;

  int sp_time_intervals_;
  int time_intervals_for_time_measurements_;

  int sp_fermionic_frequencies_;
  int hts_bosonic_frequencies_;
  int four_point_fermionic_frequencies_;

  double min_real_frequency_;
  double max_real_frequency_;
  int real_frequencies_;
  double imaginary_damping_;
};

template <typename Concurrency>
int DomainsParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(cluster_);
  buffer_size += concurrency.get_buffer_size(sp_host_);
  buffer_size += concurrency.get_buffer_size(tp_host_);
  buffer_size += concurrency.get_buffer_size(sp_time_intervals_);
  buffer_size += concurrency.get_buffer_size(time_intervals_for_time_measurements_);
  buffer_size += concurrency.get_buffer_size(sp_fermionic_frequencies_);
  buffer_size += concurrency.get_buffer_size(hts_bosonic_frequencies_);
  buffer_size += concurrency.get_buffer_size(four_point_fermionic_frequencies_);
  buffer_size += concurrency.get_buffer_size(min_real_frequency_);
  buffer_size += concurrency.get_buffer_size(max_real_frequency_);
  buffer_size += concurrency.get_buffer_size(real_frequencies_);
  buffer_size += concurrency.get_buffer_size(imaginary_damping_);

  return buffer_size;
}

template <typename Concurrency>
void DomainsParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                             int& position) const {
  concurrency.pack(buffer, buffer_size, position, cluster_);
  concurrency.pack(buffer, buffer_size, position, sp_host_);
  concurrency.pack(buffer, buffer_size, position, tp_host_);
  concurrency.pack(buffer, buffer_size, position, sp_time_intervals_);
  concurrency.pack(buffer, buffer_size, position, time_intervals_for_time_measurements_);
  concurrency.pack(buffer, buffer_size, position, sp_fermionic_frequencies_);
  concurrency.pack(buffer, buffer_size, position, hts_bosonic_frequencies_);
  concurrency.pack(buffer, buffer_size, position, four_point_fermionic_frequencies_);
  concurrency.pack(buffer, buffer_size, position, min_real_frequency_);
  concurrency.pack(buffer, buffer_size, position, max_real_frequency_);
  concurrency.pack(buffer, buffer_size, position, real_frequencies_);
  concurrency.pack(buffer, buffer_size, position, imaginary_damping_);
}

template <typename Concurrency>
void DomainsParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                               int& position) {
  concurrency.unpack(buffer, buffer_size, position, cluster_);
  concurrency.unpack(buffer, buffer_size, position, sp_host_);
  concurrency.unpack(buffer, buffer_size, position, tp_host_);
  concurrency.unpack(buffer, buffer_size, position, sp_time_intervals_);
  concurrency.unpack(buffer, buffer_size, position, time_intervals_for_time_measurements_);
  concurrency.unpack(buffer, buffer_size, position, sp_fermionic_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, hts_bosonic_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, four_point_fermionic_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, min_real_frequency_);
  concurrency.unpack(buffer, buffer_size, position, max_real_frequency_);
  concurrency.unpack(buffer, buffer_size, position, real_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, imaginary_damping_);
}
template <typename ReaderOrWriter>
void DomainsParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("domains");

    try {
      reader_or_writer.open_group("real-space-grids");

      try {
        reader_or_writer.execute("cluster", cluster_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("sp-host", sp_host_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("tp-host", tp_host_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.open_group("imaginary-time");

      try {
        reader_or_writer.execute("sp-time-intervals", sp_time_intervals_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("time-intervals-for-time-measurements",
                                 time_intervals_for_time_measurements_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.open_group("imaginary-frequency");

      try {
        reader_or_writer.execute("sp-fermionic-frequencies", sp_fermionic_frequencies_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("HTS-bosonic-frequencies", hts_bosonic_frequencies_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("four-point-fermionic-frequencies",
                                 four_point_fermionic_frequencies_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.open_group("real-frequency");

      try {
        reader_or_writer.execute("min", min_real_frequency_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("max", max_real_frequency_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("frequencies", real_frequencies_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("imaginary-damping", imaginary_damping_);
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

#endif  // DCA_PHYS_PARAMETERS_DOMAINS_PARAMETERS_HPP
