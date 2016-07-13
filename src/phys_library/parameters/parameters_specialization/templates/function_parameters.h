// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class contains parameters that define physics domains. These domains are then used to define
// 'physical' functions.

#ifndef PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_FUNCTION_PARAMETERS_H
#define PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_FUNCTION_PARAMETERS_H

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>

class function_parameters {
public:
  function_parameters(int lattice_dimension);

  /******************************************
   ***        CONCURRENCY                 ***
   ******************************************/

  template <class concurrency_type>
  int get_buffer_size(concurrency_type& concurrency);

  template <class concurrency_type>
  void pack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template <class concurrency_type>
  void unpack(concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  /******************************************
   ***        READ/WRITE                  ***
   ******************************************/

  template <class read_write_type>
  void read_write(read_write_type& read_write_obj);

  /******************************************
   ***        DATA                        ***
   ******************************************/
  double get_min_real_frequency();
  double get_max_real_frequency();
  int get_number_of_real_frequencies();
  double get_real_frequencies_off_set();

  std::vector<int> get_H_k_grid_size();

  int get_sp_time_intervals();

  int get_sp_fermionic_frequencies();
  int get_sp_bosonic_frequencies();

  std::vector<std::vector<int>> get_sp_cluster();

  int get_tp_time_intervals();

  int get_tp_fermionic_frequencies();
  int get_tp_bosonic_frequencies();

  std::vector<std::vector<int>> get_tp_cluster();

private:
  std::vector<int> H_k_grid_size;

  int sp_time_intervals;

  int sp_fermionic_frequencies;
  int sp_bosonic_frequencies;

  std::vector<std::vector<int>> sp_cluster;

  double lower_bound;
  double upper_bound;

  int nr_intervals;
  double real_axis_off_set;

  int tp_time_intervals;

  int tp_fermionic_frequencies;
  int tp_bosonic_frequencies;

  std::vector<std::vector<int>> tp_cluster;
};

function_parameters::function_parameters(int lattice_dimension)
    : H_k_grid_size(lattice_dimension, 8),

      sp_time_intervals(128),

      sp_fermionic_frequencies(256),
      sp_bosonic_frequencies(32),

      sp_cluster(0),

      lower_bound(-10.),
      upper_bound(10),

      nr_intervals(128),
      real_axis_off_set(0.01),

      tp_time_intervals(0),

      tp_fermionic_frequencies(0),
      tp_bosonic_frequencies(0),

      tp_cluster(0) {}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template <class concurrency_type>
int function_parameters::get_buffer_size(concurrency_type& concurrency) {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(H_k_grid_size);

  buffer_size += concurrency.get_buffer_size(sp_time_intervals);

  buffer_size += concurrency.get_buffer_size(sp_fermionic_frequencies);
  buffer_size += concurrency.get_buffer_size(sp_bosonic_frequencies);

  buffer_size += concurrency.get_buffer_size(sp_cluster);

  buffer_size += concurrency.get_buffer_size(lower_bound);
  buffer_size += concurrency.get_buffer_size(upper_bound);

  buffer_size += concurrency.get_buffer_size(nr_intervals);
  buffer_size += concurrency.get_buffer_size(real_axis_off_set);

  buffer_size += concurrency.get_buffer_size(tp_time_intervals);

  buffer_size += concurrency.get_buffer_size(tp_fermionic_frequencies);
  buffer_size += concurrency.get_buffer_size(tp_bosonic_frequencies);

  buffer_size += concurrency.get_buffer_size(tp_cluster);

  return buffer_size;
}

template <class concurrency_type>
void function_parameters::pack(concurrency_type& concurrency, int* buffer, int buffer_size,
                               int& position) {
  concurrency.pack(buffer, buffer_size, position, H_k_grid_size);

  concurrency.pack(buffer, buffer_size, position, sp_time_intervals);

  concurrency.pack(buffer, buffer_size, position, sp_fermionic_frequencies);
  concurrency.pack(buffer, buffer_size, position, sp_bosonic_frequencies);

  concurrency.pack(buffer, buffer_size, position, sp_cluster);

  concurrency.pack(buffer, buffer_size, position, lower_bound);
  concurrency.pack(buffer, buffer_size, position, upper_bound);
  concurrency.pack(buffer, buffer_size, position, nr_intervals);
  concurrency.pack(buffer, buffer_size, position, real_axis_off_set);

  concurrency.pack(buffer, buffer_size, position, tp_time_intervals);

  concurrency.pack(buffer, buffer_size, position, tp_fermionic_frequencies);
  concurrency.pack(buffer, buffer_size, position, tp_bosonic_frequencies);

  concurrency.pack(buffer, buffer_size, position, tp_cluster);
}

template <class concurrency_type>
void function_parameters::unpack(concurrency_type& concurrency, int* buffer, int buffer_size,
                                 int& position) {
  concurrency.unpack(buffer, buffer_size, position, H_k_grid_size);

  concurrency.unpack(buffer, buffer_size, position, sp_time_intervals);

  concurrency.unpack(buffer, buffer_size, position, sp_fermionic_frequencies);
  concurrency.unpack(buffer, buffer_size, position, sp_bosonic_frequencies);

  concurrency.unpack(buffer, buffer_size, position, sp_cluster);

  {
    concurrency.unpack(buffer, buffer_size, position, lower_bound);
    concurrency.unpack(buffer, buffer_size, position, upper_bound);
    concurrency.unpack(buffer, buffer_size, position, nr_intervals);
    concurrency.unpack(buffer, buffer_size, position, real_axis_off_set);
  }

  {
    concurrency.unpack(buffer, buffer_size, position, tp_time_intervals);

    concurrency.unpack(buffer, buffer_size, position, tp_fermionic_frequencies);
    concurrency.unpack(buffer, buffer_size, position, tp_bosonic_frequencies);

    concurrency.unpack(buffer, buffer_size, position, tp_cluster);
  }
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template <class read_write_type>
void function_parameters::read_write(read_write_type& read_write_obj) {
  try {
    read_write_obj.open_group("function-parameters");

    {
      read_write_obj.open_group("single-particle-functions");

      try {
        read_write_obj.execute("H(k) grid-size", H_k_grid_size);
      }
      catch (const std::exception& r_e) {
        std::cout << "\n not read : H(k) grid-size \n";
      }

      try {
        read_write_obj.execute("time-intervals", sp_time_intervals);
      }
      catch (const std::exception& r_e) {
      }

      try {
        read_write_obj.execute("fermionic-frequencies", sp_fermionic_frequencies);
      }
      catch (const std::exception& r_e) {
      }
      try {
        read_write_obj.execute("bosonic-frequencies", sp_bosonic_frequencies);
      }
      catch (const std::exception& r_e) {
      }

      try {
        read_write_obj.execute("sp-cluster", sp_cluster);
      }
      catch (const std::exception& r_e) {
      }

      read_write_obj.close_group();
    }

    {
      read_write_obj.open_group("two-particle-functions");

      try {
        read_write_obj.execute("time-intervals", tp_time_intervals);
      }
      catch (const std::exception& r_e) {
      }

      try {
        read_write_obj.execute("fermionic-frequencies", tp_fermionic_frequencies);
      }
      catch (const std::exception& r_e) {
      }
      try {
        read_write_obj.execute("bosonic-frequencies", tp_bosonic_frequencies);
      }
      catch (const std::exception& r_e) {
      }

      try {
        read_write_obj.execute("tp-cluster", tp_cluster);
      }
      catch (const std::exception& r_e) {
      }

      read_write_obj.close_group();
    }

    {
      read_write_obj.open_group("real-axis-functions");

      try {
        read_write_obj.execute("lower-bound", lower_bound);
      }
      catch (const std::exception& r_e) {
      }
      try {
        read_write_obj.execute("upper-bound", upper_bound);
      }
      catch (const std::exception& r_e) {
      }

      try {
        read_write_obj.execute("nr-intervals", nr_intervals);
      }
      catch (const std::exception& r_e) {
      }

      try {
        read_write_obj.execute("real-axis-off-set", real_axis_off_set);
      }
      catch (const std::exception& r_e) {
      }

      read_write_obj.close_group();
    }

    read_write_obj.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\n\t MCI-parameters defined !!  \n\n";
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

double function_parameters::get_min_real_frequency() {
  return lower_bound;
}

double function_parameters::get_max_real_frequency() {
  return upper_bound;
}

int function_parameters::get_number_of_real_frequencies() {
  return nr_intervals;
}

double function_parameters::get_real_frequencies_off_set() {
  assert(real_axis_off_set > 0);
  return real_axis_off_set;
}

std::vector<int> function_parameters::get_H_k_grid_size() {
  return H_k_grid_size;
}

int function_parameters::get_sp_time_intervals() {
  return sp_time_intervals;
}

int function_parameters::get_sp_fermionic_frequencies() {
  return sp_fermionic_frequencies;
}

int function_parameters::get_sp_bosonic_frequencies() {
  return sp_bosonic_frequencies;
}

std::vector<std::vector<int>> function_parameters::get_sp_cluster() {
  return sp_cluster;
}

int function_parameters::get_tp_time_intervals() {
  return tp_time_intervals;
}

int function_parameters::get_tp_fermionic_frequencies() {
  return tp_fermionic_frequencies;
}

int function_parameters::get_tp_bosonic_frequencies() {
  return tp_bosonic_frequencies;
}

std::vector<std::vector<int>> function_parameters::get_tp_cluster() {
  return tp_cluster;
}

#endif  // PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_FUNCTION_PARAMETERS_H
