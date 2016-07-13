// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class contains all parameters for the double-counting correction.

#ifndef PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_DOUBLE_COUNTING_PARAMETERS_H
#define PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_DOUBLE_COUNTING_PARAMETERS_H

#include <stdexcept>
#include <string>

class double_counting_parameters {
public:
  double_counting_parameters();

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

  std::string get_double_counting_method();
  double get_double_counting_correction();

  double& get_double_counting_correction_for_this_iteration();

private:
  std::string double_counting_method;
  double double_counting_correction;

  double double_counting_correction_for_this_iteration;
};

double_counting_parameters::double_counting_parameters()
    : double_counting_method("none"), double_counting_correction(0) {}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template <class concurrency_type>
int double_counting_parameters::get_buffer_size(concurrency_type& concurrency) {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(double_counting_method);
  buffer_size += concurrency.get_buffer_size(double_counting_correction);

  return buffer_size;
}

template <class concurrency_type>
void double_counting_parameters::pack(concurrency_type& concurrency, int* buffer, int buffer_size,
                                      int& position) {
  concurrency.pack(buffer, buffer_size, position, double_counting_method);
  concurrency.pack(buffer, buffer_size, position, double_counting_correction);
}

template <class concurrency_type>
void double_counting_parameters::unpack(concurrency_type& concurrency, int* buffer, int buffer_size,
                                        int& position) {
  concurrency.unpack(buffer, buffer_size, position, double_counting_method);
  concurrency.unpack(buffer, buffer_size, position, double_counting_correction);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template <class read_write_type>
void double_counting_parameters::read_write(read_write_type& read_write_obj) {
  try {
    read_write_obj.open_group("double-counting-parameters");

    try {
      read_write_obj.execute("double-counting-method", double_counting_method);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("double-counting-correction", double_counting_correction);
    }
    catch (const std::exception& r_e) {
    }

    read_write_obj.close_group();
  }
  catch (const std::exception& r_e) {
  }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

std::string double_counting_parameters::get_double_counting_method() {
  return double_counting_method;
}

double double_counting_parameters::get_double_counting_correction() {
  return double_counting_correction;
}

double& double_counting_parameters::get_double_counting_correction_for_this_iteration() {
  return double_counting_correction_for_this_iteration;
}

#endif  // PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_DOUBLE_COUNTING_PARAMETERS_H
