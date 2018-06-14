// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements the high temperature series expansion solver.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_HIGH_TEMPERATURE_SERIES_EXPANSION_SOLVER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_HIGH_TEMPERATURE_SERIES_EXPANSION_SOLVER_HPP

#include <iostream>
#include <stdexcept>
#include <string>

#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_writer.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_step/cluster_solver/high_temperature_series_expansion/series_expansion_sigma.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
class HighTemperatureSeriesExpansionSolver {
public:
  HighTemperatureSeriesExpansionSolver(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  void initialize();
  void initialize(int dca_iteration);

  void execute();

  void finalize();

  template <typename dca_info_struct_t>
  void finalize(dca_info_struct_t& dca_info_struct);

  void read(std::string filename);

  void write(std::string filename);

private:
  template <typename Writer>
  void write(Writer& writer);

private:
  parameters_type& parameters;
  MOMS_type& MOMS;

  htseries::SeriesExpansionSigma<parameters_type, MOMS_type> series_exp_obj;

  int DCA_it;
};

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
HighTemperatureSeriesExpansionSolver<device_t, parameters_type, MOMS_type>::HighTemperatureSeriesExpansionSolver(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref)
    : parameters(parameters_ref),
      MOMS(MOMS_ref),

      series_exp_obj(parameters, MOMS) {}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void HighTemperatureSeriesExpansionSolver<device_t, parameters_type, MOMS_type>::initialize() {}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void HighTemperatureSeriesExpansionSolver<device_t, parameters_type, MOMS_type>::initialize(
    int dca_iteration) {
  DCA_it = dca_iteration;
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void HighTemperatureSeriesExpansionSolver<device_t, parameters_type, MOMS_type>::execute() {
  series_exp_obj.execute(false);

  // for(int i=0; i<MOMS.Sigma.size(); ++i)
  // MOMS.Sigma(i) = series_exp_obj.get_Sigma()(i);

  for (int i = 0; i < MOMS.Sigma_lattice.size(); ++i)
    MOMS.Sigma_lattice(i) = series_exp_obj.get_Sigma()(i);

  MOMS.compute_Sigma_bands();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void HighTemperatureSeriesExpansionSolver<device_t, parameters_type, MOMS_type>::finalize() {}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
template <typename dca_info_struct_t>
void HighTemperatureSeriesExpansionSolver<device_t, parameters_type, MOMS_type>::finalize(
    dca_info_struct_t& /*dca_info_struct*/) {
  finalize();
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
void HighTemperatureSeriesExpansionSolver<device_t, parameters_type, MOMS_type>::write(
    std::string file_name) {
  std::cout << "\n\n\t\t start writing " << file_name << "\n\n";

  const std::string& output_format = parameters.get_output_format();

  if (output_format == "JSON") {
    dca::io::JSONWriter writer;
    writer.open_file(file_name);

    parameters.write(writer);
    MOMS.write(writer);
    this->write(writer);

    writer.close_file();
  }

  else if (output_format == "HDF5") {
    dca::io::HDF5Writer writer;
    writer.open_file(file_name);

    parameters.write(writer);
    MOMS.write(writer);
    this->write(writer);

    writer.close_file();
  }

  else
    throw std::logic_error(__FUNCTION__);
}

template <dca::linalg::DeviceType device_t, class parameters_type, class MOMS_type>
template <typename Writer>
void HighTemperatureSeriesExpansionSolver<device_t, parameters_type, MOMS_type>::write(
    Writer& /*writer*/) {
  // writer.open_group("functions");

  // series_exp_obj.write(writer);

  // writer.close_group();
}

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_HIGH_TEMPERATURE_SERIES_EXPANSION_HIGH_TEMPERATURE_SERIES_EXPANSION_SOLVER_HPP
