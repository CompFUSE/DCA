// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_MODEL_SPECILIZATIONS_MODEL_TIGHT_BINDING_PARAMETERS_H
#define PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_MODEL_SPECILIZATIONS_MODEL_TIGHT_BINDING_PARAMETERS_H

#include "phys_library/parameters/parameters_specialization/templates/model_parameters.h"

#include <iostream>
#include <stdexcept>

#include "phys_library/parameters/models/tight_binding_model.h"

template <typename lattice_type>
class model_parameters<tight_binding_model<lattice_type>> {
public:
  model_parameters();

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

  double get_t();
  double get_t_prime();

  double get_U();

  double get_V();
  double get_V_prime();

private:
  double t;
  double t_prime;

  double U;

  double V;
  double V_prime;
};

template <typename lattice_type>
model_parameters<tight_binding_model<lattice_type>>::model_parameters()
    : t(1.),
      t_prime(0.),

      U(4.),

      V(0.),
      V_prime(0.) {}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template <typename lattice_type>
template <class concurrency_type>
int model_parameters<tight_binding_model<lattice_type>>::get_buffer_size(concurrency_type& concurrency) {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t);
  buffer_size += concurrency.get_buffer_size(t_prime);

  buffer_size += concurrency.get_buffer_size(U);

  buffer_size += concurrency.get_buffer_size(V);
  buffer_size += concurrency.get_buffer_size(V_prime);

  return buffer_size;
}

template <typename lattice_type>
template <class concurrency_type>
void model_parameters<tight_binding_model<lattice_type>>::pack(concurrency_type& concurrency,
                                                               int* buffer, int buffer_size,
                                                               int& position) {
  concurrency.pack(buffer, buffer_size, position, t);
  concurrency.pack(buffer, buffer_size, position, t_prime);

  concurrency.pack(buffer, buffer_size, position, U);

  concurrency.pack(buffer, buffer_size, position, V);
  concurrency.pack(buffer, buffer_size, position, V_prime);
}

template <typename lattice_type>
template <class concurrency_type>
void model_parameters<tight_binding_model<lattice_type>>::unpack(concurrency_type& concurrency,
                                                                 int* buffer, int buffer_size,
                                                                 int& position) {
  concurrency.unpack(buffer, buffer_size, position, t);
  concurrency.unpack(buffer, buffer_size, position, t_prime);

  concurrency.unpack(buffer, buffer_size, position, U);

  concurrency.unpack(buffer, buffer_size, position, V);
  concurrency.unpack(buffer, buffer_size, position, V_prime);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template <typename lattice_type>
template <class read_write_type>
void model_parameters<tight_binding_model<lattice_type>>::read_write(read_write_type& read_write_obj) {
  try {
    read_write_obj.open_group("2D-Hubbard-model");

    try {
      read_write_obj.execute("t", t);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("t-prime", t_prime);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("U", U);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("V", V);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("V-prime", V_prime);
    }
    catch (const std::exception& r_e) {
    }

    read_write_obj.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\n\t 2D-Hubbard-model parameters defined !!\n\n";
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

template <typename lattice_type>
double model_parameters<tight_binding_model<lattice_type>>::get_t() {
  return t;
}

template <typename lattice_type>
double model_parameters<tight_binding_model<lattice_type>>::get_t_prime() {
  return t_prime;
}

template <typename lattice_type>
double model_parameters<tight_binding_model<lattice_type>>::get_U() {
  return U;
}

template <typename lattice_type>
double model_parameters<tight_binding_model<lattice_type>>::get_V() {
  return V;
}

template <typename lattice_type>
double model_parameters<tight_binding_model<lattice_type>>::get_V_prime() {
  return V_prime;
}

#endif  // PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_MODEL_SPECILIZATIONS_MODEL_TIGHT_BINDING_PARAMETERS_H
