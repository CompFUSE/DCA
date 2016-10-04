// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the high temperature series expansion solver.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_HIGH_TEMPERATURE_SERIES_EXPANSION_SOLVER_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_HIGH_TEMPERATURE_SERIES_EXPANSION_SOLVER_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_template.h"

#include <iostream>
#include <stdexcept>
#include <string>

#include "comp_library/IO_library/IO.hpp"
#include "dca/linalg/device_type.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_series_expansion/series_expansion_sigma.h"

namespace DCA {

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
class cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type> {
public:
  cluster_solver(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  void initialize();
  void initialize(int dca_iteration);

  void execute();

  void finalize();

  template <typename dca_info_struct_t>
  void finalize(dca_info_struct_t& dca_info_struct);

  void read(std::string filename);

  void write(std::string filename);

private:
  template <IO::FORMAT DATA_FORMAT>
  void write(IO::writer<DATA_FORMAT>& writer);

private:
  parameters_type& parameters;
  MOMS_type& MOMS;

  SERIES_EXPANSION::series_expansion<parameters_type, MOMS_type> series_exp_obj;

  int DCA_it;
};

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::cluster_solver(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref)
    : parameters(parameters_ref),
      MOMS(MOMS_ref),

      series_exp_obj(parameters, MOMS) {}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::initialize() {}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::initialize(
    int dca_iteration) {
  DCA_it = dca_iteration;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::execute() {
  series_exp_obj.execute(false);

  // for(int i=0; i<MOMS.Sigma.size(); ++i)
  // MOMS.Sigma(i) = series_exp_obj.get_Sigma()(i);

  for (int i = 0; i < MOMS.Sigma_lattice.size(); ++i)
    MOMS.Sigma_lattice(i) = series_exp_obj.get_Sigma()(i);

  MOMS.compute_Sigma_bands();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::finalize() {}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
template <typename dca_info_struct_t>
void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::finalize(
    dca_info_struct_t& /*dca_info_struct*/) {
  finalize();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::write(
    std::string file_name) {
  IO::FORMAT FORMAT = parameters.get_output_format();

  std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

  switch (FORMAT) {
    case IO::JSON: {
      IO::writer<IO::JSON> writer;
      {
        writer.open_file(file_name);

        parameters.write(writer);
        MOMS.write(writer);
        this->write(writer);

        writer.close_file();
      }
    } break;

    case IO::HDF5: {
      IO::writer<IO::HDF5> writer;
      {
        writer.open_file(file_name);

        parameters.write(writer);
        MOMS.write(writer);
        this->write(writer);

        writer.close_file();
      }
    } break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
template <IO::FORMAT DATA_FORMAT>
void cluster_solver<HIGH_TEMPERATURE_SERIES, device_t, parameters_type, MOMS_type>::write(
    IO::writer<DATA_FORMAT>& /*writer*/) {
  // writer.open_group("functions");

  // series_exp_obj.write(writer);

  // writer.close_group();
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_HIGH_TEMPERATURE_SERIES_EXPANSION_SOLVER_H
