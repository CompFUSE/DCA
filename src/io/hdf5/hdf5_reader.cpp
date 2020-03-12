// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements hdf5_reader.hpp.

#include "dca/io/hdf5/hdf5_reader.hpp"

#include <fstream>
#include <stdexcept>

namespace dca {
namespace io {
// dca::io::

HDF5Reader::~HDF5Reader() {
  if (file_)
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

  file_ = std::make_unique<H5::H5File>(file_name.c_str(), H5F_ACC_RDONLY);
}

void HDF5Reader::close_file() {
  file_->close();
  file_.release();
}

std::string HDF5Reader::get_path() {
  std::string path = "/";

  for (size_t i = 0; i < paths_.size(); i++) {
    path = path + paths_[i];

    if (i < paths_.size() - 1)
      path = path + "/";
  }

  return path;
}

bool HDF5Reader::execute(const std::string& name, std::string& value) {
  std::string full_name = get_path() + "/" + name;

  if (!exists(full_name)) {
    return false;
  }

  auto dims = readSize(full_name);
  assert(dims.size() == 1);
  value.resize(dims.at(0));

  read(full_name, HDF5_TYPE<char>::get_PredType(), value.data());
  return true;
}

bool HDF5Reader::execute(const std::string& name, std::vector<std::string>& value) {
  std::string full_name = get_path() + "/" + name;
  if (!exists(full_name)) {
    return false;
  }

  open_group(name);

  int size = -1;
  execute("size", size);

  value.resize(size);

  open_group("data");

  for (size_t l = 0; l < value.size(); l++) {
    open_group(std::to_string(l));

    execute(std::to_string(l), value[l]);

    close_group();
  }

  close_group();
  close_group();

  return true;
}

void HDF5Reader::read(const std::string& name, H5::PredType type, void* data) const {
  H5::DataSet dataset = file_->openDataSet(name.c_str());
  dataset.read(data, type);
}

bool HDF5Reader::exists(const std::string& name) const {
  auto code = H5Gget_objinfo(file_->getId(), name.c_str(), 0, NULL);
  return code == 0;
}

std::vector<hsize_t> HDF5Reader::readSize(const std::string& name) const {
  H5::DataSet dataset = file_->openDataSet(name.c_str());
  H5::DataSpace dataspace = dataset.getSpace();

  int n_dims = dataspace.getSimpleExtentNdims();
  std::vector<hsize_t> dims(n_dims);
  dataspace.getSimpleExtentDims(dims.data(), nullptr);

  return dims;
}

}  // namespace io
}  // namespace dca
