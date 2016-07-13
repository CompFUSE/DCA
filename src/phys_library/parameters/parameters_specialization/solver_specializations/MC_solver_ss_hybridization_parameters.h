// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class spezializes the MC solver parameters for the CT-AUX QMC method.

#ifndef PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_SOLVER_SPECILIZATIONS_MC_SOLVER_SS_HYBRIDIZATION_PARAMETERS_H
#define PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_SOLVER_SPECILIZATIONS_MC_SOLVER_SS_HYBRIDIZATION_PARAMETERS_H

#include "phys_library/parameters/parameters_specialization/templates/MC_solver_parameters.h"

#include <iostream>
#include <stdexcept>

template <>
class MC_solver_parameters<DCA::SS_CT_HYB> {
public:
  MC_solver_parameters();

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

  int get_Sigma_tail_cutoff();

  double get_steps_per_sweep();
  double get_swaps_per_sweep();
  double get_shifts_per_sweep();

private:
  int Sigma_tail_cutoff;

  double steps_per_sweep;
  double swaps_per_sweep;
  double shifts_per_sweep;
};

MC_solver_parameters<DCA::SS_CT_HYB>::MC_solver_parameters()
    :

      Sigma_tail_cutoff(0),

      steps_per_sweep(0.4),
      swaps_per_sweep(0.1),
      shifts_per_sweep(0.4) {}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template <class concurrency_type>
int MC_solver_parameters<DCA::SS_CT_HYB>::get_buffer_size(concurrency_type& concurrency) {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(Sigma_tail_cutoff);

  buffer_size += concurrency.get_buffer_size(steps_per_sweep);
  buffer_size += concurrency.get_buffer_size(swaps_per_sweep);
  buffer_size += concurrency.get_buffer_size(shifts_per_sweep);

  return buffer_size;
}

template <class concurrency_type>
void MC_solver_parameters<DCA::SS_CT_HYB>::pack(concurrency_type& concurrency, int* buffer,
                                                int buffer_size, int& position) {
  concurrency.pack(buffer, buffer_size, position, Sigma_tail_cutoff);

  concurrency.pack(buffer, buffer_size, position, steps_per_sweep);
  concurrency.pack(buffer, buffer_size, position, swaps_per_sweep);
  concurrency.pack(buffer, buffer_size, position, shifts_per_sweep);
}

template <class concurrency_type>
void MC_solver_parameters<DCA::SS_CT_HYB>::unpack(concurrency_type& concurrency, int* buffer,
                                                  int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, Sigma_tail_cutoff);

  concurrency.unpack(buffer, buffer_size, position, steps_per_sweep);
  concurrency.unpack(buffer, buffer_size, position, swaps_per_sweep);
  concurrency.unpack(buffer, buffer_size, position, shifts_per_sweep);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template <class read_write_type>
void MC_solver_parameters<DCA::SS_CT_HYB>::read_write(read_write_type& read_write_obj) {
  try {
    read_write_obj.open_group("SS-CT-HYB-solver");

    try {
      read_write_obj.execute("Sigma-tail-cutoff", Sigma_tail_cutoff);
    }
    catch (const std::exception& r_e) {
    }

    try {
      read_write_obj.execute("steps-per-sweep", steps_per_sweep);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("swaps-per-sweep", swaps_per_sweep);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("shifts-per-sweep", shifts_per_sweep);
    }
    catch (const std::exception& r_e) {
    }

    read_write_obj.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\n\t SS-CT-HYB-solver-parameters defined !!  \n\n";
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

/******************************************
 ***        DATA                        ***
 ******************************************/

int MC_solver_parameters<DCA::SS_CT_HYB>::get_Sigma_tail_cutoff() {
  return Sigma_tail_cutoff;
}

double MC_solver_parameters<DCA::SS_CT_HYB>::get_steps_per_sweep() {
  return steps_per_sweep;
}

double MC_solver_parameters<DCA::SS_CT_HYB>::get_swaps_per_sweep() {
  return swaps_per_sweep;
}

double MC_solver_parameters<DCA::SS_CT_HYB>::get_shifts_per_sweep() {
  return shifts_per_sweep;
}

#endif  // PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_SOLVER_SPECILIZATIONS_MC_SOLVER_SS_HYBRIDIZATION_PARAMETERS_H
