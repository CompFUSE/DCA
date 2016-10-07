// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the filename parameters.

#ifndef DCA_PHYS_PARAMETERS_FILENAME_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_FILENAME_PARAMETERS_HPP

#include <iostream>
#include <stdexcept>
#include <string>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class FilenameParameters {
public:
  FilenameParameters()
      : directory_("./"),
        output_format_("JSON"),
        output_("output.json"),
        profiling_("prof_data.txt"),
        spectrum_("data_spectrum.json"),
        susceptibilities_("data_susceptibilities.json"),
        vertex_("vertex_filename.py"),
        ED_output_("output_ED.json"),
        CPE_output_("output_CPE.json"),
        QMC_output_("output_QMC.json"),
        dump_lattice_Self_energy_("false"),
        dump_cluster_Greens_functions_("false") {}

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
  std::string get_file_extension() const {
    if (output_format_ == "HDF5")
      return ".hdf5";
    else if (output_format_ == "JSON")
      return ".json";
    else
      throw std::logic_error(__FUNCTION__);
  }
  const std::string& get_output_file_name() const {
    return output_;
  }
  const std::string& get_profiling_file_name() const {
    return profiling_;
  }
  const std::string& get_spectrum_file_name() const {
    return spectrum_;
  }
  const std::string& get_susceptibilities_file_name() const {
    return susceptibilities_;
  }
  const std::string& get_vertex_file_name() const {
    return vertex_;
  }
  const std::string& get_ED_output_file_name() const {
    return ED_output_;
  }
  const std::string& get_CPE_output_file_name() const {
    return CPE_output_;
  }
  const std::string& get_QMC_output_file_name() const {
    return QMC_output_;
  }
  bool dump_lattice_Self_energy() const {
    return (dump_lattice_Self_energy_ == "true");
  }
  bool dump_cluster_Greens_functions() const {
    return (dump_cluster_Greens_functions_ == "true");
  }

private:
  std::string directory_;
  std::string output_format_;
  std::string output_;
  std::string profiling_;
  std::string spectrum_;
  std::string susceptibilities_;
  std::string vertex_;
  std::string ED_output_;
  std::string CPE_output_;
  std::string QMC_output_;
  std::string dump_lattice_Self_energy_;
  std::string dump_cluster_Greens_functions_;
};

template <typename Concurrency>
int FilenameParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(directory_);
  buffer_size += concurrency.get_buffer_size(output_format_);
  buffer_size += concurrency.get_buffer_size(output_);
  buffer_size += concurrency.get_buffer_size(profiling_);
  buffer_size += concurrency.get_buffer_size(spectrum_);
  buffer_size += concurrency.get_buffer_size(susceptibilities_);
  buffer_size += concurrency.get_buffer_size(vertex_);
  buffer_size += concurrency.get_buffer_size(ED_output_);
  buffer_size += concurrency.get_buffer_size(CPE_output_);
  buffer_size += concurrency.get_buffer_size(QMC_output_);
  buffer_size += concurrency.get_buffer_size(dump_lattice_Self_energy_);
  buffer_size += concurrency.get_buffer_size(dump_cluster_Greens_functions_);

  return buffer_size;
}

template <typename Concurrency>
void FilenameParameters::pack(const Concurrency& concurrency, int* buffer, int buffer_size,
                              int& position) const {
  concurrency.pack(buffer, buffer_size, position, directory_);
  concurrency.pack(buffer, buffer_size, position, output_format_);
  concurrency.pack(buffer, buffer_size, position, output_);
  concurrency.pack(buffer, buffer_size, position, profiling_);
  concurrency.pack(buffer, buffer_size, position, spectrum_);
  concurrency.pack(buffer, buffer_size, position, susceptibilities_);
  concurrency.pack(buffer, buffer_size, position, vertex_);
  concurrency.pack(buffer, buffer_size, position, ED_output_);
  concurrency.pack(buffer, buffer_size, position, CPE_output_);
  concurrency.pack(buffer, buffer_size, position, QMC_output_);
  concurrency.pack(buffer, buffer_size, position, dump_lattice_Self_energy_);
  concurrency.pack(buffer, buffer_size, position, dump_cluster_Greens_functions_);
}

template <typename Concurrency>
void FilenameParameters::unpack(const Concurrency& concurrency, int* buffer, int buffer_size,
                                int& position) {
  concurrency.unpack(buffer, buffer_size, position, directory_);
  concurrency.unpack(buffer, buffer_size, position, output_format_);
  concurrency.unpack(buffer, buffer_size, position, output_);
  concurrency.unpack(buffer, buffer_size, position, profiling_);
  concurrency.unpack(buffer, buffer_size, position, spectrum_);
  concurrency.unpack(buffer, buffer_size, position, susceptibilities_);
  concurrency.unpack(buffer, buffer_size, position, vertex_);
  concurrency.unpack(buffer, buffer_size, position, ED_output_);
  concurrency.unpack(buffer, buffer_size, position, CPE_output_);
  concurrency.unpack(buffer, buffer_size, position, QMC_output_);
  concurrency.unpack(buffer, buffer_size, position, dump_lattice_Self_energy_);
  concurrency.unpack(buffer, buffer_size, position, dump_cluster_Greens_functions_);
}

template <typename ReaderOrWriter>
void FilenameParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("filename-parameters");

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
      reader_or_writer.execute("output-file", output_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("profiling-file", profiling_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("spectrum-file", spectrum_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("susceptibilities-file", susceptibilities_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("plot-vertex-file", vertex_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("output-ED", ED_output_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("output-CPE", CPE_output_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("output-QMC", QMC_output_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("dump-lattice-self-energy", dump_lattice_Self_energy_);
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
    std::cout << "\nNo filename parameters defined!" << std::endl;
    ;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_FILENAME_PARAMETERS_HPP
