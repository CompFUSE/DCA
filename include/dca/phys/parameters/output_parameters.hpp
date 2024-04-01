// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This class reads, stores, and writes the output parameters.

#ifndef DCA_PHYS_PARAMETERS_OUTPUT_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_OUTPUT_PARAMETERS_HPP

#include <stdexcept>
#include <string>
#include "dca/io/io_types.hpp"

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class OutputParameters {
public:
  OutputParameters()
      : directory_("./"),
	autoresume_(false),
#ifdef DCA_HAVE_ADIOS2
        output_format_("HDF5"),
        g4_output_format_(""),
        filename_dca_("dca.hdf5"),
        filename_analysis_("sofqomega.bp"),
#else
        output_format_("HDF5"),
        filename_dca_("dca.hdf5"),
        filename_analysis_("sofqomega.bp"),
#endif
        directory_config_read_(""),
        directory_config_write_(""),
        filename_ed_("ed.hdf5"),
        filename_qmc_("qmc.hdf5"),
        filename_profiling_("profiling.json"),
#ifdef DCA_HAVE_ADIOS2
        autoresume_filename_{filename_dca_},
#else
        autoresume_filename_{filename_dca_ + std::string{".tmp"}},
#endif
        dump_lattice_self_energy_(false),
        dump_cluster_Greens_functions_(false),
        dump_Gamma_lattice_(false),
        dump_chi_0_lattice_(false) {
  }

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const std::string& get_directory() const {
    return directory_;
  }

  bool autoresume() const {
    return autoresume_;
  }
  const std::string& get_autoresume_filename() const {
    return autoresume_filename_;
  }
  const std::string& get_output_format() const {
    return output_format_;
  }
  const std::string& get_g4_output_format() const {
    return g4_output_format_;
  }
  const std::string& get_directory_config_read() const {
    return directory_config_read_;
  }
  const std::string& get_directory_config_write() const {
    return directory_config_write_;
  }
  const std::string& get_filename_g4() const {
    return filename_g4_;
  }
  const std::string& get_filename_dca() const {
    return filename_dca_;
  }
  const std::string& get_filename_analysis() const {
    return filename_analysis_;
  }
  std::string getAppropriateFilenameAnalysis(int rank = 0) const {
    std::size_t extension_start = filename_analysis_.rfind('.');
    std::string cleaned_filename{filename_analysis_.substr(0, extension_start) + std::to_string(rank)};
    cleaned_filename += extensionFromIOType(io::stringToIOType(get_output_format()));
    return cleaned_filename;
  }
  const std::string& get_filename_ed() const {
    return filename_ed_;
  }
  const std::string& get_filename_qmc() const {
    return filename_qmc_;
  }
  const std::string& get_filename_profiling() const {
    return filename_profiling_;
  }
  bool dump_lattice_self_energy() const {
    return dump_lattice_self_energy_;
  }
  bool dump_cluster_Greens_functions() const {
    return dump_cluster_Greens_functions_;
  }
  bool dump_Gamma_lattice() const {
    return dump_Gamma_lattice_;
  }
  bool dump_chi_0_lattice() const {
    return dump_chi_0_lattice_;
  }
  bool dump_every_iteration() const {
    return dump_every_iteration_;
  }

private:
  std::string directory_;
  bool autoresume_; // = false;
  std::string output_format_;
  std::string g4_output_format_;
  std::string filename_g4_;
  std::string filename_dca_;
  std::string filename_analysis_;
  std::string directory_config_read_;
  std::string directory_config_write_;
  std::string filename_ed_;
  std::string filename_qmc_;
  std::string filename_profiling_;
  std::string autoresume_filename_;
  bool dump_lattice_self_energy_;
  bool dump_cluster_Greens_functions_;
  bool dump_Gamma_lattice_;
  bool dump_chi_0_lattice_;
  bool dump_every_iteration_ = false;
};

template <typename Concurrency>
int OutputParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(directory_);
  buffer_size += concurrency.get_buffer_size(autoresume_);
  buffer_size += concurrency.get_buffer_size(output_format_);
  buffer_size += concurrency.get_buffer_size(g4_output_format_);
  buffer_size += concurrency.get_buffer_size(filename_g4_);
  buffer_size += concurrency.get_buffer_size(filename_dca_);
  buffer_size += concurrency.get_buffer_size(filename_analysis_);
  buffer_size += concurrency.get_buffer_size(directory_config_read_);
  buffer_size += concurrency.get_buffer_size(directory_config_write_);
  buffer_size += concurrency.get_buffer_size(filename_ed_);
  buffer_size += concurrency.get_buffer_size(filename_qmc_);
  buffer_size += concurrency.get_buffer_size(filename_profiling_);
  buffer_size += concurrency.get_buffer_size(autoresume_filename_);
  buffer_size += concurrency.get_buffer_size(dump_lattice_self_energy_);
  buffer_size += concurrency.get_buffer_size(dump_cluster_Greens_functions_);
  buffer_size += concurrency.get_buffer_size(dump_Gamma_lattice_);
  buffer_size += concurrency.get_buffer_size(dump_chi_0_lattice_);
  buffer_size += concurrency.get_buffer_size(dump_every_iteration_);

  return buffer_size;
}

template <typename Concurrency>
void OutputParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                            int& position) const {
  concurrency.pack(buffer, buffer_size, position, directory_);
  concurrency.pack(buffer, buffer_size, position, autoresume_);
  concurrency.pack(buffer, buffer_size, position, output_format_);
  concurrency.pack(buffer, buffer_size, position, g4_output_format_);
  concurrency.pack(buffer, buffer_size, position, filename_g4_);
  concurrency.pack(buffer, buffer_size, position, filename_dca_);
  concurrency.pack(buffer, buffer_size, position, filename_analysis_);
  concurrency.pack(buffer, buffer_size, position, directory_config_read_);
  concurrency.pack(buffer, buffer_size, position, directory_config_write_);
  concurrency.pack(buffer, buffer_size, position, filename_ed_);
  concurrency.pack(buffer, buffer_size, position, filename_qmc_);
  concurrency.pack(buffer, buffer_size, position, filename_profiling_);
  concurrency.pack(buffer, buffer_size, position, autoresume_filename_);
  concurrency.pack(buffer, buffer_size, position, dump_lattice_self_energy_);
  concurrency.pack(buffer, buffer_size, position, dump_cluster_Greens_functions_);
  concurrency.pack(buffer, buffer_size, position, dump_Gamma_lattice_);
  concurrency.pack(buffer, buffer_size, position, dump_chi_0_lattice_);
  concurrency.pack(buffer, buffer_size, position, dump_every_iteration_);
}

template <typename Concurrency>
void OutputParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                              int& position) {
  concurrency.unpack(buffer, buffer_size, position, directory_);
  concurrency.unpack(buffer, buffer_size, position, autoresume_);
  concurrency.unpack(buffer, buffer_size, position, output_format_);
  concurrency.unpack(buffer, buffer_size, position, g4_output_format_);
  concurrency.unpack(buffer, buffer_size, position, filename_g4_);
  concurrency.unpack(buffer, buffer_size, position, filename_dca_);
  concurrency.unpack(buffer, buffer_size, position, filename_analysis_);
  concurrency.unpack(buffer, buffer_size, position, directory_config_read_);
  concurrency.unpack(buffer, buffer_size, position, directory_config_write_);
  concurrency.unpack(buffer, buffer_size, position, filename_ed_);
  concurrency.unpack(buffer, buffer_size, position, filename_qmc_);
  concurrency.unpack(buffer, buffer_size, position, filename_profiling_);
  concurrency.unpack(buffer, buffer_size, position, autoresume_filename_);
  concurrency.unpack(buffer, buffer_size, position, dump_lattice_self_energy_);
  concurrency.unpack(buffer, buffer_size, position, dump_cluster_Greens_functions_);
  concurrency.unpack(buffer, buffer_size, position, dump_Gamma_lattice_);
  concurrency.unpack(buffer, buffer_size, position, dump_chi_0_lattice_);
  concurrency.unpack(buffer, buffer_size, position, dump_every_iteration_);
}

template <typename ReaderOrWriter>
void OutputParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  auto try_to_read_or_write = [&](const std::string& name, auto& obj) {
    try {
      reader_or_writer.execute(name, obj);
    }
    catch (const std::exception& r_e) {
    }
  };

  try {
    reader_or_writer.open_group("output");

    try_to_read_or_write("directory", directory_);
    try_to_read_or_write("autoresume", autoresume_);
    try_to_read_or_write("output-format", output_format_);
    try_to_read_or_write("g4-output-format", g4_output_format_);
    try_to_read_or_write("filename-g4", filename_g4_);
    try_to_read_or_write("filename-dca", filename_dca_);
    try_to_read_or_write("filename-analysis", filename_analysis_);
    try_to_read_or_write("directory-config-read", directory_config_read_);
    try_to_read_or_write("directory-config-write", directory_config_write_);
    try_to_read_or_write("filename-ed", filename_ed_);
    try_to_read_or_write("filename-qmc", filename_qmc_);
    try_to_read_or_write("filename-profiling", filename_profiling_);
    try_to_read_or_write("autoresume-filename", autoresume_filename_);
    try_to_read_or_write("dump-lattice-self-energy", dump_lattice_self_energy_);
    try_to_read_or_write("dump-cluster-Greens-functions", dump_cluster_Greens_functions_);
    try_to_read_or_write("dump-Gamma-lattice", dump_Gamma_lattice_);
    try_to_read_or_write("dump-chi-0-lattice", dump_chi_0_lattice_);
    try_to_read_or_write("dump-every-iteration", dump_every_iteration_);

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }

  // Check here for known invalid input, if we know options conflict don't just crash while running.

  if (dump_every_iteration_) {
    io ::IOType io_type = io::stringToIOType(output_format_);
    switch (io_type) {
      // google test cannot handle these cases being combined.?
      case io::IOType::ADIOS2:
	break;
      case io::IOType::HDF5:
        break;
      case io::IOType::JSON:
        throw std::runtime_error("JSON output format does not support dump-ever-iteration.");
    }
  }
}

}  // namespace params
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_PARAMETERS_OUTPUT_PARAMETERS_HPP
