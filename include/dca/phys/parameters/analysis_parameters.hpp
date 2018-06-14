// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class reads, stores, and writes the analysis parameters.

#ifndef DCA_PHYS_PARAMETERS_ANALYSIS_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_ANALYSIS_PARAMETERS_HPP

#include <iostream>
#include <stdexcept>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class AnalysisParameters {
public:
  AnalysisParameters()
      : symmetrize_Gamma_(true),
        Gamma_deconvolution_cut_off_(0.5),
        project_onto_crystal_harmonics_(false),
        projection_cut_off_radius_(1.5) {}

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

private:
  bool symmetrize_Gamma_;
  double Gamma_deconvolution_cut_off_;
  bool project_onto_crystal_harmonics_;
  double projection_cut_off_radius_;
};

template <typename Concurrency>
int AnalysisParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(symmetrize_Gamma_);
  buffer_size += concurrency.get_buffer_size(Gamma_deconvolution_cut_off_);
  buffer_size += concurrency.get_buffer_size(project_onto_crystal_harmonics_);
  buffer_size += concurrency.get_buffer_size(projection_cut_off_radius_);

  return buffer_size;
}

template <typename Concurrency>
void AnalysisParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                              int& position) const {
  concurrency.pack(buffer, buffer_size, position, symmetrize_Gamma_);
  concurrency.pack(buffer, buffer_size, position, Gamma_deconvolution_cut_off_);
  concurrency.pack(buffer, buffer_size, position, project_onto_crystal_harmonics_);
  concurrency.pack(buffer, buffer_size, position, projection_cut_off_radius_);
}

template <typename Concurrency>
void AnalysisParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                                int& position) {
  concurrency.unpack(buffer, buffer_size, position, symmetrize_Gamma_);
  concurrency.unpack(buffer, buffer_size, position, Gamma_deconvolution_cut_off_);
  concurrency.unpack(buffer, buffer_size, position, project_onto_crystal_harmonics_);
  concurrency.unpack(buffer, buffer_size, position, projection_cut_off_radius_);
}
template <typename ReaderOrWriter>
void AnalysisParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("analysis");

    try {
      reader_or_writer.execute("symmetrize-Gamma", symmetrize_Gamma_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("Gamma-deconvolution-cut-off", Gamma_deconvolution_cut_off_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("project-onto-crystal-harmonics", project_onto_crystal_harmonics_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("projection-cut-off-radius", projection_cut_off_radius_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_ANALYSIS_PARAMETERS_HPP
