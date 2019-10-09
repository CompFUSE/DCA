// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements hdf5_reader.hpp.

#include "dca/io/hdf5/hdf5_reader.hpp"

#include <fstream>
#include <stdexcept>

namespace dca {
namespace io {
// dca::io::

HDF5Reader::~HDF5Reader() {
  if (my_file != NULL)
    close_file();
}

void HDF5Reader::open_file(std::string file_name) {
  {  // check whether the file exists ...
    std::wifstream tmp(file_name.c_str());

    if (!tmp or !tmp.good() or tmp.bad()) {
      throw std::runtime_error("Cannot open file : " + file_name);
    }
    else if (verbose_) {
      std::cout << "\n\n\topening file : " << file_name << "\n";
    }
  }

  my_file = new H5::H5File(file_name.c_str(), H5F_ACC_RDONLY);
}

void HDF5Reader::close_file() {
  delete my_file;
  my_file = NULL;
}

std::string HDF5Reader::get_path() {
  std::string path = "/";

  for (size_t i = 0; i < my_paths.size(); i++) {
    path = path + my_paths[i];

    if (i < my_paths.size() - 1)
      path = path + "/";
  }

  return path;
}

bool HDF5Reader::execute(std::string name,
                         std::string& value)  //, H5File& file, std::string path)
{
  std::string full_name = get_path() + "/" + name;

  try {
    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    value.resize(dataset.getInMemDataSize(), 'a');

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<char>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &value[0]);
    return true;
  }
  catch (const H5::FileIException& err) {
    std::cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    return false;
  }
}

bool HDF5Reader::execute(std::string name, std::vector<std::string>& value) {
  open_group(name);
  open_group("data");
  bool success = true;

  try {
    int size = -1;
    execute("size", size);

    value.resize(size);

    for (size_t l = 0; l < value.size(); l++) {
      std::stringstream ss;
      ss << get_path() << "/" << l;

      H5::DataSet dataset = my_file->openDataSet(ss.str().c_str());

      value[l].resize(dataset.getInMemDataSize(), 'a');

      H5::DataSpace dataspace = dataset.getSpace();

      H5Dread(dataset.getId(), HDF5_TYPE<char>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
              &value[l][0]);
    }
  }
  catch (const H5::FileIException& err) {
    std::cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    success = false;
  }

  close_group();
  close_group();
  return success;
}

}  // namespace io
}  // namespace dca
