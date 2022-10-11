// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// This class reads, stores, and writes the analysis parameters.

#ifndef DCA_PHYS_PARAMETERS_ANALYSIS_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_ANALYSIS_PARAMETERS_HPP

#include <iostream>
#include <stdexcept>
#include <vector>
#include "dca/phys/four_point_type.hpp"

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class AnalysisParameters {
public:
  AnalysisParameters(int dimension)
      : symmetrize_Gamma_(true),
        Gamma_deconvolution_cut_off_(0.5),
        project_onto_crystal_harmonics_(false),
        projection_cut_off_radius_(1.5),
	q_host_(dimension, std::vector<int>(dimension,0))
		{
    for (int i = 0; i < dimension; ++i) {
      q_host_[i][i] = 1;
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

  bool symmetrize_Gamma() const {
    return symmetrize_Gamma_;
  }
  double get_Gamma_deconvolution_cut_off() const {
    return Gamma_deconvolution_cut_off_;
  }
  bool project_onto_crystal_harmonics() const {
    return project_onto_crystal_harmonics_;
  }
  double get_projection_cut_off_radius() const {
    return projection_cut_off_radius_;
  }
  const std::vector<std::vector<int>>& get_q_host() const {
    return q_host_;
  }

private:
  bool symmetrize_Gamma_ = true;
  double Gamma_deconvolution_cut_off_ = 0.5;
  bool project_onto_crystal_harmonics_ = false;
  double projection_cut_off_radius_ = 1.5;
  FourPointType g4_channel_ = FourPointType::PARTICLE_HOLE_MAGNETIC;
  std::vector<std::vector<int>> q_host_;
};

template <typename Concurrency>
int AnalysisParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(symmetrize_Gamma_);
  buffer_size += concurrency.get_buffer_size(Gamma_deconvolution_cut_off_);
  buffer_size += concurrency.get_buffer_size(project_onto_crystal_harmonics_);
  buffer_size += concurrency.get_buffer_size(projection_cut_off_radius_);
  buffer_size += concurrency.get_buffer_size(g4_channel_);
  buffer_size += concurrency.get_buffer_size(q_host_);

  return buffer_size;
}

template <typename Concurrency>
void AnalysisParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                              int& position) const {
  concurrency.pack(buffer, buffer_size, position, symmetrize_Gamma_);
  concurrency.pack(buffer, buffer_size, position, Gamma_deconvolution_cut_off_);
  concurrency.pack(buffer, buffer_size, position, project_onto_crystal_harmonics_);
  concurrency.pack(buffer, buffer_size, position, projection_cut_off_radius_);
  concurrency.pack(buffer, buffer_size, position, g4_channel_);
  concurrency.pack(buffer, buffer_size, position, q_host_);
}

template <typename Concurrency>
void AnalysisParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                                int& position) {
  concurrency.unpack(buffer, buffer_size, position, symmetrize_Gamma_);
  concurrency.unpack(buffer, buffer_size, position, Gamma_deconvolution_cut_off_);
  concurrency.unpack(buffer, buffer_size, position, project_onto_crystal_harmonics_);
  concurrency.unpack(buffer, buffer_size, position, projection_cut_off_radius_);
  concurrency.unpack(buffer, buffer_size, position, g4_channel_);
  concurrency.unpack(buffer, buffer_size, position, q_host_);
}

template <typename ReaderOrWriter>
void AnalysisParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  if (reader_or_writer.open_group("analysis")) {
    reader_or_writer.execute("symmetrize-Gamma", symmetrize_Gamma_);
    reader_or_writer.execute("Gamma-deconvolution-cut-off", Gamma_deconvolution_cut_off_);
    reader_or_writer.execute("project-onto-crystal-harmonics", project_onto_crystal_harmonics_);
    reader_or_writer.execute("projection-cut-off-radius", projection_cut_off_radius_);
    if (reader_or_writer.is_reader) {
      std::string g4_channel_name;
      reader_or_writer.execute("g4-channel", g4_channel_name);
      g4_channel_ = stringToFourPointType(g4_channel_name);
    }
    else {
      std::string g4_channel_name = toString(g4_channel_);
      reader_or_writer.execute("g4-channel", g4_channel_name);
    }
    try {
      reader_or_writer.execute("q-host", q_host_);
    }
    catch (const std::exception& r_e) {
    }
    reader_or_writer.close_group();
  }
}

}  // namespace params
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_PARAMETERS_ANALYSIS_PARAMETERS_HPP
