// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the function parameters.

#ifndef DCA_PHYS_PARAMETERS_FUNCTION_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_FUNCTION_PARAMETERS_HPP

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class FunctionParameters {
public:
  FunctionParameters(int lattice_dimension)
      : H_k_grid_size_(lattice_dimension, 8),

        sp_time_intervals_(128),
        sp_fermionic_frequencies_(256),
        sp_bosonic_frequencies_(32),
        sp_cluster_(0),

        lower_bound_(-10.),
        upper_bound_(10),
        nr_intervals_(128),
        real_axis_off_set_(0.01),

        tp_time_intervals_(0),
        tp_fermionic_frequencies_(0),
        tp_bosonic_frequencies_(0),
        tp_cluster_(0) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const std::vector<int>& get_H_k_grid_size() const {
    return H_k_grid_size_;
  }

  int get_sp_time_intervals() const {
    return sp_time_intervals_;
  }
  int get_sp_fermionic_frequencies() const {
    return sp_fermionic_frequencies_;
  }
  int get_sp_bosonic_frequencies() const {
    return sp_bosonic_frequencies_;
  }
  const std::vector<std::vector<int>>& get_sp_cluster() const {
    return sp_cluster_;
  }

  double get_min_real_frequency() const {
    return lower_bound_;
  }
  double get_max_real_frequency() const {
    return upper_bound_;
  }
  int get_number_of_real_frequencies() const {
    return nr_intervals_;
  }
  double get_real_frequencies_off_set() const {
    assert(real_axis_off_set_ > 0);
    return real_axis_off_set_;
  }

  int get_tp_time_intervals() const {
    return tp_time_intervals_;
  }
  int get_tp_fermionic_frequencies() const {
    return tp_fermionic_frequencies_;
  }
  int get_tp_bosonic_frequencies() const {
    return tp_bosonic_frequencies_;
  }
  const std::vector<std::vector<int>>& get_tp_cluster() const {
    return tp_cluster_;
  }

private:
  std::vector<int> H_k_grid_size_;

  int sp_time_intervals_;
  int sp_fermionic_frequencies_;
  int sp_bosonic_frequencies_;
  std::vector<std::vector<int>> sp_cluster_;

  double lower_bound_;
  double upper_bound_;
  int nr_intervals_;
  double real_axis_off_set_;

  int tp_time_intervals_;
  int tp_fermionic_frequencies_;
  int tp_bosonic_frequencies_;
  std::vector<std::vector<int>> tp_cluster_;
};

template <typename Concurrency>
int FunctionParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(H_k_grid_size_);

  buffer_size += concurrency.get_buffer_size(sp_time_intervals_);
  buffer_size += concurrency.get_buffer_size(sp_fermionic_frequencies_);
  buffer_size += concurrency.get_buffer_size(sp_bosonic_frequencies_);
  buffer_size += concurrency.get_buffer_size(sp_cluster_);

  buffer_size += concurrency.get_buffer_size(lower_bound_);
  buffer_size += concurrency.get_buffer_size(upper_bound_);
  buffer_size += concurrency.get_buffer_size(nr_intervals_);
  buffer_size += concurrency.get_buffer_size(real_axis_off_set_);

  buffer_size += concurrency.get_buffer_size(tp_time_intervals_);
  buffer_size += concurrency.get_buffer_size(tp_fermionic_frequencies_);
  buffer_size += concurrency.get_buffer_size(tp_bosonic_frequencies_);
  buffer_size += concurrency.get_buffer_size(tp_cluster_);

  return buffer_size;
}

template <typename Concurrency>
void FunctionParameters::pack(const Concurrency& concurrency, int* buffer, int buffer_size,
                              int& position) const {
  concurrency.pack(buffer, buffer_size, position, H_k_grid_size_);

  concurrency.pack(buffer, buffer_size, position, sp_time_intervals_);
  concurrency.pack(buffer, buffer_size, position, sp_fermionic_frequencies_);
  concurrency.pack(buffer, buffer_size, position, sp_bosonic_frequencies_);
  concurrency.pack(buffer, buffer_size, position, sp_cluster_);

  concurrency.pack(buffer, buffer_size, position, lower_bound_);
  concurrency.pack(buffer, buffer_size, position, upper_bound_);
  concurrency.pack(buffer, buffer_size, position, nr_intervals_);
  concurrency.pack(buffer, buffer_size, position, real_axis_off_set_);

  concurrency.pack(buffer, buffer_size, position, tp_time_intervals_);
  concurrency.pack(buffer, buffer_size, position, tp_fermionic_frequencies_);
  concurrency.pack(buffer, buffer_size, position, tp_bosonic_frequencies_);
  concurrency.pack(buffer, buffer_size, position, tp_cluster_);
}

template <typename Concurrency>
void FunctionParameters::unpack(const Concurrency& concurrency, int* buffer, int buffer_size,
                                int& position) {
  concurrency.unpack(buffer, buffer_size, position, H_k_grid_size_);

  concurrency.unpack(buffer, buffer_size, position, sp_time_intervals_);
  concurrency.unpack(buffer, buffer_size, position, sp_fermionic_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, sp_bosonic_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, sp_cluster_);

  concurrency.unpack(buffer, buffer_size, position, lower_bound_);
  concurrency.unpack(buffer, buffer_size, position, upper_bound_);
  concurrency.unpack(buffer, buffer_size, position, nr_intervals_);
  concurrency.unpack(buffer, buffer_size, position, real_axis_off_set_);

  concurrency.unpack(buffer, buffer_size, position, tp_time_intervals_);
  concurrency.unpack(buffer, buffer_size, position, tp_fermionic_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, tp_bosonic_frequencies_);
  concurrency.unpack(buffer, buffer_size, position, tp_cluster_);
}

template <typename ReaderOrWriter>
void FunctionParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("function-parameters");

    {
      reader_or_writer.open_group("single-particle-functions");

      try {
        reader_or_writer.execute("H(k) grid-size", H_k_grid_size_);
      }
      catch (const std::exception& r_e) {
        std::cout << "\nnot read: H(k) grid-size\n" << std::endl;
      }

      try {
        reader_or_writer.execute("time-intervals", sp_time_intervals_);
      }
      catch (const std::exception& r_e) {
      }

      try {
        reader_or_writer.execute("fermionic-frequencies", sp_fermionic_frequencies_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("bosonic-frequencies", sp_bosonic_frequencies_);
      }
      catch (const std::exception& r_e) {
      }

      try {
        reader_or_writer.execute("sp-cluster", sp_cluster_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }

    {
      reader_or_writer.open_group("two-particle-functions");

      try {
        reader_or_writer.execute("time-intervals", tp_time_intervals_);
      }
      catch (const std::exception& r_e) {
      }

      try {
        reader_or_writer.execute("fermionic-frequencies", tp_fermionic_frequencies_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("bosonic-frequencies", tp_bosonic_frequencies_);
      }
      catch (const std::exception& r_e) {
      }

      try {
        reader_or_writer.execute("tp-cluster", tp_cluster_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }

    {
      reader_or_writer.open_group("real-axis-functions");

      try {
        reader_or_writer.execute("lower-bound", lower_bound_);
      }
      catch (const std::exception& r_e) {
      }
      try {
        reader_or_writer.execute("upper-bound", upper_bound_);
      }
      catch (const std::exception& r_e) {
      }

      try {
        reader_or_writer.execute("nr-intervals", nr_intervals_);
      }
      catch (const std::exception& r_e) {
      }

      try {
        reader_or_writer.execute("real-axis-off-set", real_axis_off_set_);
      }
      catch (const std::exception& r_e) {
      }

      reader_or_writer.close_group();
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\nNo function parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

}  // params
}  // phys
}  // dca
#endif  // DCA_PHYS_PARAMETERS_FUNCTIONPARAMETERS_HPP
