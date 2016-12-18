// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class reads, stores, and writes the output parameters.

#ifndef DCA_PHYS_PARAMETERS_OUTPUT_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_OUTPUT_PARAMETERS_HPP

#include <stdexcept>
#include <string>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class OutputParameters {
public:
  OutputParameters()
      : directory_("./"),
        output_format_("HDF5"),
        filename_dca_("dca.hdf5"),
        filename_analysis_("analysis.hdf5"),
        filename_ed_("ed.hdf5"),
        filename_qmc_("qmc.hdf5"),
        filename_profiling_("profiling.json"),
        dump_lattice_self_energy_(false),
        dump_cluster_Greens_functions_(false) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const std::string& get_directory() const {
    return directory_;
  }
  const std::string& get_output_format() const {
    return output_format_;
  }
  const std::string& get_filename_dca() const {
    return filename_dca_;
  }
  const std::string& get_filename_analysis() const {
    return filename_analysis_;
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

private:
  std::string directory_;
  std::string output_format_;
  std::string filename_dca_;
  std::string filename_analysis_;
  std::string filename_ed_;
  std::string filename_qmc_;
  std::string filename_profiling_;
  bool dump_lattice_self_energy_;
  bool dump_cluster_Greens_functions_;
};

template <typename Concurrency>
int OutputParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(directory_);
  buffer_size += concurrency.get_buffer_size(output_format_);
  buffer_size += concurrency.get_buffer_size(filename_dca_);
  buffer_size += concurrency.get_buffer_size(filename_analysis_);
  buffer_size += concurrency.get_buffer_size(filename_ed_);
  buffer_size += concurrency.get_buffer_size(filename_qmc_);
  buffer_size += concurrency.get_buffer_size(filename_profiling_);
  buffer_size += concurrency.get_buffer_size(dump_lattice_self_energy_);
  buffer_size += concurrency.get_buffer_size(dump_cluster_Greens_functions_);

  return buffer_size;
}

template <typename Concurrency>
void OutputParameters::pack(const Concurrency& concurrency, int* buffer, int buffer_size,
                            int& position) const {
  concurrency.pack(buffer, buffer_size, position, directory_);
  concurrency.pack(buffer, buffer_size, position, output_format_);
  concurrency.pack(buffer, buffer_size, position, filename_dca_);
  concurrency.pack(buffer, buffer_size, position, filename_analysis_);
  concurrency.pack(buffer, buffer_size, position, filename_ed_);
  concurrency.pack(buffer, buffer_size, position, filename_qmc_);
  concurrency.pack(buffer, buffer_size, position, filename_profiling_);
  concurrency.pack(buffer, buffer_size, position, dump_lattice_self_energy_);
  concurrency.pack(buffer, buffer_size, position, dump_cluster_Greens_functions_);
}

template <typename Concurrency>
void OutputParameters::unpack(const Concurrency& concurrency, int* buffer, int buffer_size,
                              int& position) {
  concurrency.unpack(buffer, buffer_size, position, directory_);
  concurrency.unpack(buffer, buffer_size, position, output_format_);
  concurrency.unpack(buffer, buffer_size, position, filename_dca_);
  concurrency.unpack(buffer, buffer_size, position, filename_analysis_);
  concurrency.unpack(buffer, buffer_size, position, filename_ed_);
  concurrency.unpack(buffer, buffer_size, position, filename_qmc_);
  concurrency.unpack(buffer, buffer_size, position, filename_profiling_);
  concurrency.unpack(buffer, buffer_size, position, dump_lattice_self_energy_);
  concurrency.unpack(buffer, buffer_size, position, dump_cluster_Greens_functions_);
}

template <typename ReaderOrWriter>
void OutputParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("output");

    try {
      reader_or_writer.execute("directory", directory_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("output-format", output_format_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("filename-dca", filename_dca_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("filename-analysis", filename_analysis_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("filename-ed", filename_ed_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("filename-qmc", filename_qmc_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("filename-profiling", filename_profiling_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("dump-lattice-self-energy", dump_lattice_self_energy_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("dump-cluster-Greens-function", dump_cluster_Greens_functions_);
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

#endif  // DCA_PHYS_PARAMETERS_OUTPUT_PARAMETERS_HPP
