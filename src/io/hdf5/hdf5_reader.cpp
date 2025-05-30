// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// This file implements hdf5_reader.hpp.

#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/filesystem.hpp"
#include "ModernStringUtils.hpp"
#include "hdf5.h"
#include <fstream>
#include <stdexcept>

namespace dca {
namespace io {
// dca::io::

HDF5Reader::~HDF5Reader() {
  if (file_)
    close_file();
}

long HDF5Reader::getStepCount() {
  long steps;
  bool has_steps = execute("steps", steps);
  if (!has_steps) {
    is_legacy_ = true;
    std::cerr << "Legacy DCA hdf5 with no step data read.\n";
    return -1;
  }
  return steps;
}

void HDF5Reader::open_file(std::string file_name) {
  if (file_)
    throw std::logic_error(__FUNCTION__);

  try {
    if (filesystem::exists(file_name))
      file_ = std::make_unique<H5::H5File>(file_name.c_str(), H5F_ACC_RDONLY);
    else
      throw std::runtime_error("Cannot open file: " + file_name + " to read, it is not found!");

    if (verbose_)
      std::cout << "\n\n\topening file : " << file_name << "\n";
  }
  catch (const std::exception& ex) {
    throw std::runtime_error("Cannot open file : " + file_name);
  }
}

void HDF5Reader::close_file() {
  file_->close();
  file_.reset();
  in_step_ = false;
  step_ = 0;
  paths_.clear();
}

std::string HDF5Reader::get_path() const {
  std::string path = "/";

  for (size_t i = 0; i < paths_.size(); i++) {
    path = path + paths_[i];

    if (i < paths_.size() - 1)
      path = path + "/";
  }

  return path;
}

void HDF5Reader::begin_step() {
  if (is_legacy_)
    return;
  if (in_step_)
    throw std::runtime_error("HDF5Reader::begin_step() called while already in step!");
  in_step_ = true;
  std::string step_group{"step_" + std::to_string(++step_ - 1)};
  paths_.push_back(step_group);
}

void HDF5Reader::end_step() {
  if (is_legacy_)
    return;
  if (!in_step_)
    throw std::runtime_error("HDF5Reader::end_step() called while not in step!");
  paths_.clear();
  in_step_ = false;
}

bool HDF5Reader::execute(const std::string& name, std::string& value) {
  std::string full_name;
  if (!buildCheckedFullName(name, full_name))
    return false;

  H5::DataSet dataset = file_->openDataSet(full_name.c_str());
  const auto type = dataset.getDataType();

  const auto size = type.getSize();
  value.resize(size);

  dataset.read(value.data(), type);

  // Null string case.
  if (value == std::string{0})
    value = "";

  return true;
}

bool HDF5Reader::execute(const std::string& name, std::vector<std::string>& value) {
  std::string full_name = get_path();
  if (full_name.size() < 1)
    full_name += name;
  else
    full_name += "/" + name;
  if (!exists(full_name)) {
    return false;
  }

  H5::DataSet dataset = file_->openDataSet(name.c_str());
  auto size = readSize(full_name)[0];
  auto s_type = H5::StrType(H5::PredType::C_S1, H5T_VARIABLE);

  std::vector<char*> data(size);
  dataset.read(data.data(), s_type);

  value.resize(size);
  for (int i = 0; i < size; ++i) {
    value[i] = data[i];
  }

  // clean memory
  dataset.vlenReclaim(data.data(), s_type, dataset.getSpace());

  return true;
}

void HDF5Reader::read(const std::string& name, H5::DataType type, void* data) const {
  H5::DataSet dataset = file_->openDataSet(name.c_str());
  dataset.read(data, type);
}

bool HDF5Reader::exists(const std::string& name) const {
  // can't supress exception from C++ library if more than one group level is checked at once.
  return file_->nameExists(name);
}

bool HDF5Reader::buildCheckedFullName(const std::string& name, std::string& full_name) const {
  full_name = get_path();
  auto parts = modernstrutil::split(name, "/");
  std::string part_wise_name;
  bool is_found = false;
  for (auto part : parts) {
    part_wise_name += std::string{"/"} + std::string{part};
    std::string tmp_string = full_name + part_wise_name;
    is_found = exists(tmp_string);
    if (!is_found)
      return false;
  }
  full_name += part_wise_name;
  return true;
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
