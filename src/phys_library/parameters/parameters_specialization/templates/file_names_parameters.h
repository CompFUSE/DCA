// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class contains all parameters for reading and writing (input/output) files.

#ifndef PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_FILE_NAMES_PARAMETERS_H
#define PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_FILE_NAMES_PARAMETERS_H

#include <iostream>
#include <stdexcept>
#include <string>

#include "comp_library/IO_library/IO_types.h"

class file_names_parameters {
public:
  file_names_parameters();

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

  std::string& get_directory();

  std::string& get_input_file_name();
  std::string& get_output_file_name();

  std::string& get_profiling_file_name();

  std::string& get_spectrum_file_name();
  std::string& get_susceptibilities_file_name();

  std::string& get_vertex_file_name();

  std::string& get_ED_output_file_name();
  std::string& get_CPE_output_file_name();
  std::string& get_QMC_output_file_name();

  IO::FORMAT get_output_format();
  std::string get_file_extension();

  bool dump_lattice_Self_energy();
  bool dump_cluster_Greens_functions();

private:
  std::string directory_filename;

  std::string input_filename;
  std::string output_filename;

  std::string profiling_filename;

  std::string spectrum_filename;
  std::string susceptibilities_filename;

  std::string vertex_filename;

  std::string ED_output_file;
  std::string CPE_output_file;
  std::string QMC_output_file;

  std::string format_str;

  std::string dump_lattice_Self_energy_str;
  std::string dump_cluster_Greens_functions_str;
};

file_names_parameters::file_names_parameters()
    : directory_filename("./"),

      input_filename("input.json"),
      output_filename("output.json"),

      profiling_filename("prof_data.txt"),

      spectrum_filename("data_spectrum.json"),
      susceptibilities_filename("data_susceptibilities.json"),

      vertex_filename("vertex_filename.py"),

      ED_output_file("output_ED.json"),
      CPE_output_file("output_CPE.json"),
      QMC_output_file("output_QMC.json"),

      format_str("JSON"),

      dump_lattice_Self_energy_str("false"),
      dump_cluster_Greens_functions_str("false") {}

/******************************************
 ***        CONCURRENCY                 ***
 ******************************************/

template <class concurrency_type>
int file_names_parameters::get_buffer_size(concurrency_type& concurrency) {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(directory_filename);

  buffer_size += concurrency.get_buffer_size(input_filename);
  buffer_size += concurrency.get_buffer_size(output_filename);

  buffer_size += concurrency.get_buffer_size(profiling_filename);

  buffer_size += concurrency.get_buffer_size(spectrum_filename);
  buffer_size += concurrency.get_buffer_size(susceptibilities_filename);

  buffer_size += concurrency.get_buffer_size(vertex_filename);

  buffer_size += concurrency.get_buffer_size(ED_output_file);
  buffer_size += concurrency.get_buffer_size(CPE_output_file);
  buffer_size += concurrency.get_buffer_size(QMC_output_file);

  buffer_size += concurrency.get_buffer_size(format_str);

  buffer_size += concurrency.get_buffer_size(dump_lattice_Self_energy_str);
  buffer_size += concurrency.get_buffer_size(dump_cluster_Greens_functions_str);

  return buffer_size;
}

template <class concurrency_type>
void file_names_parameters::pack(concurrency_type& concurrency, int* buffer, int buffer_size,
                                 int& position) {
  concurrency.pack(buffer, buffer_size, position, directory_filename);

  concurrency.pack(buffer, buffer_size, position, input_filename);
  concurrency.pack(buffer, buffer_size, position, output_filename);

  concurrency.pack(buffer, buffer_size, position, profiling_filename);

  concurrency.pack(buffer, buffer_size, position, spectrum_filename);
  concurrency.pack(buffer, buffer_size, position, susceptibilities_filename);

  concurrency.pack(buffer, buffer_size, position, vertex_filename);

  concurrency.pack(buffer, buffer_size, position, ED_output_file);
  concurrency.pack(buffer, buffer_size, position, CPE_output_file);
  concurrency.pack(buffer, buffer_size, position, QMC_output_file);

  concurrency.pack(buffer, buffer_size, position, format_str);

  concurrency.pack(buffer, buffer_size, position, dump_lattice_Self_energy_str);
  concurrency.pack(buffer, buffer_size, position, dump_cluster_Greens_functions_str);
}

template <class concurrency_type>
void file_names_parameters::unpack(concurrency_type& concurrency, int* buffer, int buffer_size,
                                   int& position) {
  concurrency.unpack(buffer, buffer_size, position, directory_filename);

  concurrency.unpack(buffer, buffer_size, position, input_filename);
  concurrency.unpack(buffer, buffer_size, position, output_filename);

  concurrency.unpack(buffer, buffer_size, position, profiling_filename);

  concurrency.unpack(buffer, buffer_size, position, spectrum_filename);
  concurrency.unpack(buffer, buffer_size, position, susceptibilities_filename);

  concurrency.unpack(buffer, buffer_size, position, vertex_filename);

  concurrency.unpack(buffer, buffer_size, position, ED_output_file);
  concurrency.unpack(buffer, buffer_size, position, CPE_output_file);
  concurrency.unpack(buffer, buffer_size, position, QMC_output_file);

  concurrency.unpack(buffer, buffer_size, position, format_str);

  concurrency.unpack(buffer, buffer_size, position, dump_lattice_Self_energy_str);
  concurrency.unpack(buffer, buffer_size, position, dump_cluster_Greens_functions_str);
}

/******************************************
 ***        READ/WRITE                  ***
 ******************************************/

template <class read_write_type>
void file_names_parameters::read_write(read_write_type& read_write_obj) {
  try {
    read_write_obj.open_group("filename-parameters");

    try {
      read_write_obj.execute("directory", directory_filename);
    }
    catch (const std::exception& r_e) {
    }

    try {
      read_write_obj.execute("output-file", output_filename);
    }
    catch (const std::exception& r_e) {
    }

    try {
      read_write_obj.execute("profiling-file", profiling_filename);
    }
    catch (const std::exception& r_e) {
    }

    try {
      read_write_obj.execute("spectrum-file", spectrum_filename);
    }
    catch (const std::exception& r_e) {
    }

    try {
      read_write_obj.execute("susceptibilities-file", susceptibilities_filename);
    }
    catch (const std::exception& r_e) {
    }

    try {
      read_write_obj.execute("plot-vertex-file", vertex_filename);
    }
    catch (const std::exception& r_e) {
    }

    try {
      read_write_obj.execute("output-format", format_str);
    }
    catch (const std::exception& r_e) {
    }

    try {
      read_write_obj.execute("output-ED", ED_output_file);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("output-CPE", CPE_output_file);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("output-QMC", QMC_output_file);
    }
    catch (const std::exception& r_e) {
    }

    try {
      read_write_obj.execute("dump-lattice-self-energy", dump_lattice_Self_energy_str);
    }
    catch (const std::exception& r_e) {
    }
    try {
      read_write_obj.execute("dump-cluster-Greens-function", dump_cluster_Greens_functions_str);
    }
    catch (const std::exception& r_e) {
    }

    read_write_obj.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\n\t NO filename-parameters with output-file defined !!  \n\n";
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

std::string& file_names_parameters::get_directory() {
  return directory_filename;
}

std::string& file_names_parameters::get_input_file_name() {
  return input_filename;
}

std::string& file_names_parameters::get_output_file_name() {
  return output_filename;
}

std::string& file_names_parameters::get_profiling_file_name() {
  return profiling_filename;
}

std::string& file_names_parameters::get_spectrum_file_name() {
  return spectrum_filename;
}

std::string& file_names_parameters::get_susceptibilities_file_name() {
  return susceptibilities_filename;
}

std::string& file_names_parameters::get_vertex_file_name() {
  return vertex_filename;
}

IO::FORMAT file_names_parameters::get_output_format() {
  if (format_str == "HDF5")
    return IO::HDF5;

  return IO::JSON;
}

std::string& file_names_parameters::get_ED_output_file_name() {
  return ED_output_file;
}

std::string& file_names_parameters::get_CPE_output_file_name() {
  return CPE_output_file;
}

std::string& file_names_parameters::get_QMC_output_file_name() {
  return QMC_output_file;
}

std::string file_names_parameters::get_file_extension() {
  if (format_str == "HDF5")
    return ".hdf5";

  return ".json";
}

bool file_names_parameters::dump_lattice_Self_energy() {
  if (dump_lattice_Self_energy_str == "true")
    return true;

  return false;
}

bool file_names_parameters::dump_cluster_Greens_functions() {
  if (dump_cluster_Greens_functions_str == "true")
    return true;

  return false;
}

#endif  // PHYS_LIBRARY_PARAMETERS_PARAMETERS_SPECIALIZATION_TEMPLATES_FILE_NAMES_PARAMETERS_H
