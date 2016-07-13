// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class contains all parameters for the Monte Carlo solver.

#ifndef PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_MC_SOLVER_PARAMETERS_H
#define PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_MC_SOLVER_PARAMETERS_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_name.hpp"

template <DCA::ClusterSolverName solver_name>
class MC_solver_parameters {
public:
  MC_solver_parameters();

  /******************************************
  ***        PHYSICS                      ***
  ******************************************/
  double get_K_CT_AUX();
  int get_K_PHANI();
  int get_initial_matrix_size();

  /******************************************
   ***        CONCURRENCY                 ***
   ******************************************/

  template <class concurrency_type>
  int get_buffer_size(const concurrency_type& concurrency) const;

  template <class concurrency_type>
  void pack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  template <class concurrency_type>
  void unpack(const concurrency_type& concurrency, int* buffer, int buffer_size, int& position);

  /******************************************
   ***        READ/WRITE                  ***
   ******************************************/

  template <class stream_type>
  void to_JSON(stream_type& ss);

  template <class JSON_reader_type>
  void from_JSON(JSON_reader_type& reader);

  template <class read_write_type>
  void read_write(read_write_type& read_write_obj);
};

template <DCA::ClusterSolverName solver_name>
MC_solver_parameters<solver_name>::MC_solver_parameters() {}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template <DCA::ClusterSolverName solver_name>
template <class concurrency_type>
int MC_solver_parameters<solver_name>::get_buffer_size(const concurrency_type& /*concurrency*/) const {
  return 0;
}

template <DCA::ClusterSolverName solver_name>
template <class concurrency_type>
void MC_solver_parameters<solver_name>::pack(const concurrency_type& /*concurrency*/, int* /*buffer*/,
                                             int /*buffer_size*/, int& /*position*/) {}

template <DCA::ClusterSolverName solver_name>
template <class concurrency_type>
void MC_solver_parameters<solver_name>::unpack(const concurrency_type& /*concurrency*/,
                                               int* /*buffer*/, int /*buffer_size*/,
                                               int& /*position*/) {}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template <DCA::ClusterSolverName solver_name>
template <class stream_type>
void MC_solver_parameters<solver_name>::to_JSON(stream_type& /*ss*/) {}

template <DCA::ClusterSolverName solver_name>
template <class JSON_reader_type>
void MC_solver_parameters<solver_name>::from_JSON(JSON_reader_type& /*reader*/) {}

#endif  // PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_MC_SOLVER_PARAMETERS_H
