// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//
// Description

#ifndef PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_MODEL_SPECILIZATIONS_MODEL_2D_4BAND_MODEL_PARAMETERS_H
#define PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_MODEL_SPECILIZATIONS_MODEL_2D_4BAND_MODEL_PARAMETERS_H

#include "phys_library/parameters/parameters_specialization/templates/model_parameters.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "phys_library/parameters/models/tight_binding_model.h"
#include "phys_library/parameters/models/analytic_Hamiltonians/lattices/2D_4band_lattice.h"

template <typename dca_point_group_t>
class model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>> {
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

  double get_ei0();
  double get_eb0();
  double get_t0();

  double get_ei1();
  double get_eb1();
  double get_t1();

  double get_U0();
  double get_U1();
  double get_V();
  double get_V_prime();

private:
  double ei0;
  double eb0;
  double t0;
  double ei1;
  double eb1;
  double t1;
  double U0;
  double U1;
  double V;
  double V_prime;
};

template <typename dca_point_group_t>
model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::model_parameters()
    : ei0(0), eb0(0), t0(0), ei1(0), eb1(0), t1(0), U0(0), U1(0), V(0), V_prime(0) {}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template <typename dca_point_group_t>
template <class concurrency_type>
int model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::get_buffer_size(
    concurrency_type& concurrency) {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(ei0);
  buffer_size += concurrency.get_buffer_size(eb0);
  buffer_size += concurrency.get_buffer_size(t0);
  buffer_size += concurrency.get_buffer_size(ei1);
  buffer_size += concurrency.get_buffer_size(eb1);
  buffer_size += concurrency.get_buffer_size(t1);

  buffer_size += concurrency.get_buffer_size(U0);
  buffer_size += concurrency.get_buffer_size(U1);
  buffer_size += concurrency.get_buffer_size(V);
  buffer_size += concurrency.get_buffer_size(V_prime);

  return buffer_size;
}

template <typename dca_point_group_t>
template <class concurrency_type>
void model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::pack(
    concurrency_type& concurrency, int* buffer, int buffer_size, int& position) {
  concurrency.pack(buffer, buffer_size, position, ei0);
  concurrency.pack(buffer, buffer_size, position, eb0);
  concurrency.pack(buffer, buffer_size, position, t0);
  concurrency.pack(buffer, buffer_size, position, ei1);
  concurrency.pack(buffer, buffer_size, position, eb1);
  concurrency.pack(buffer, buffer_size, position, t1);

  concurrency.pack(buffer, buffer_size, position, U0);
  concurrency.pack(buffer, buffer_size, position, U1);
  concurrency.pack(buffer, buffer_size, position, V);
  concurrency.pack(buffer, buffer_size, position, V_prime);
}

template <typename dca_point_group_t>
template <class concurrency_type>
void model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::unpack(
    concurrency_type& concurrency, int* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, ei0);
  concurrency.unpack(buffer, buffer_size, position, eb0);
  concurrency.unpack(buffer, buffer_size, position, t0);
  concurrency.unpack(buffer, buffer_size, position, ei1);
  concurrency.unpack(buffer, buffer_size, position, eb1);
  concurrency.unpack(buffer, buffer_size, position, t1);

  concurrency.unpack(buffer, buffer_size, position, U0);
  concurrency.unpack(buffer, buffer_size, position, U1);
  concurrency.unpack(buffer, buffer_size, position, V);
  concurrency.unpack(buffer, buffer_size, position, V_prime);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template <typename dca_point_group_t>
template <class read_write_type>
void model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::read_write(
    read_write_type& read_write_obj) {
  try {
    read_write_obj.open_group("fourband-model");

    try {
      read_write_obj.execute("ei0", ei0);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("eb0", eb0);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("t0", t0);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("ei1", ei1);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("eb1", eb1);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("t1", t1);
    }
    catch (const std::exception& r_e) {
    }

    try {
      read_write_obj.execute("U0", U0);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("U1", U1);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("V", V);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("V_prime", V_prime);
    }
    catch (const std::exception& r_e) {
    }

    read_write_obj.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\n\t fourband-model parameters defined !!\n\n";
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  {
    std::stringstream ss;

    ss << "\n\n";
    ss << "\t fourband-model : \n";
    ss << "\t--------------- \n\n";
    ss << "\t\t ei0  : " << ei0 << "\n";
    ss << "\t\t eb0  : " << eb0 << "\n";
    ss << "\t\t  t0  : " << t0 << "\n";
    ss << "\t\t ei1  : " << ei1 << "\n";
    ss << "\t\t eb1  : " << eb1 << "\n";
    ss << "\t\t  t1  : " << t1 << "\n";
    ss << "\t\t  U0  : " << U0 << "\n";
    ss << "\t\t  U1  : " << U1 << "\n";
    ss << "\t\t  V   : " << V << "\n";
    ss << "\t\t  V_prime : " << V_prime << "\n";
    ss << "\n\n";

    std::cout << ss.str();
  }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

template <typename dca_point_group_t>
double model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::get_ei0() {
  return ei0;
}

template <typename dca_point_group_t>
double model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::get_eb0() {
  return eb0;
}

template <typename dca_point_group_t>
double model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::get_t0() {
  return t0;
}

template <typename dca_point_group_t>
double model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::get_ei1() {
  return ei1;
}

template <typename dca_point_group_t>
double model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::get_eb1() {
  return eb1;
}

template <typename dca_point_group_t>
double model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::get_t1() {
  return t1;
}

template <typename dca_point_group_t>
double model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::get_U0() {
  return U0;
}

template <typename dca_point_group_t>
double model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::get_U1() {
  return U1;
}

template <typename dca_point_group_t>
double model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::get_V() {
  return V;
}

template <typename dca_point_group_t>
double model_parameters<tight_binding_model<fourband_lattice<dca_point_group_t>>>::get_V_prime() {
  return V_prime;
}

#endif  // PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_MODEL_SPECILIZATIONS_MODEL_2D_4BAND_MODEL_PARAMETERS_H
