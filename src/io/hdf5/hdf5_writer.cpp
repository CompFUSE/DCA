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
// This file implements hdf5_writer.hpp.

#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/filesystem.hpp"

#include <stdexcept>

namespace dca {
namespace io {
// dca::io::

HDF5Writer::~HDF5Writer() {
  if (file_)
    close_file();
}

void HDF5Writer::open_file(std::string file_name, bool overwrite) {
  if (file_)
    throw std::logic_error(__FUNCTION__);

  if (overwrite) {
    file_ = std::make_unique<H5::H5File>(file_name.c_str(), H5F_ACC_TRUNC);
  }
  else {
    if (filesystem::exists(file_name))
      file_ = std::make_unique<H5::H5File>(file_name.c_str(), H5F_ACC_RDWR);
    else
      file_ = std::make_unique<H5::H5File>(file_name.c_str(), H5F_ACC_EXCL);
  }

  file_id_ = file_->getId();
}

void HDF5Writer::close_file() {
  if (file_) {
    // file_->flush(H5F_SCOPE_LOCAL);
    file_->close();
    file_.release();
  }
}

void HDF5Writer::open_group(std::string name) {
  my_paths_.push_back(name);
  const std::string path = get_path();

  if (!exists(path)) {
    file_->createGroup(path.c_str());
  }
}

void HDF5Writer::close_group() {
  my_paths_.pop_back();
}

std::string HDF5Writer::get_path() {
  std::string path = "/";

  for (size_t i = 0; i < my_paths_.size(); i++) {
    path += my_paths_[i];

    if (i < my_paths_.size() - 1)
      path += "/";
  }

  return path;
}

void HDF5Writer::execute(const std::string& name,
                         const std::string& value)  //, H5File& file, std::string path)
{
  // HDF5 does not allow a datatype of size 0.
  if (value.size() == 0)
    return execute(name, std::string{0});

  std::string full_name = get_path() + '/' + name;

  // String type.
  H5::StrType datatype(H5::PredType::C_S1, value.size());

  write(full_name, std::vector<hsize_t>{1}, datatype, value.data());
}

void HDF5Writer::execute(const std::string& name,
                         const std::vector<std::string>& value)  //, H5File& file, std::string path)
{
  if (value.size() == 0)
    return;

  const std::string full_name = get_path() + "/" + name;

  auto s_type = H5::StrType(H5::PredType::C_S1, H5T_VARIABLE);

  std::vector<const char*> data;
  for (const auto& s : value) {
    data.push_back(s.data());
  }

  write(full_name, std::vector<hsize_t>{data.size()}, s_type, data.data());
}

H5::DataSet HDF5Writer::write(const std::string& name, const std::vector<hsize_t>& dims,
                              H5::DataType type, const void* data) {
  if (exists(name)) {
    H5::DataSet dataset = file_->openDataSet(name.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    size_check_.resize(dims.size());
    dataspace.getSimpleExtentDims(size_check_.data(), nullptr);

    if (size_check_ != dims) {
      throw(std::length_error("Object size different than HDF5 dataset."));
    }

    dataset.write(data, type, dataspace, H5P_DEFAULT);
    return dataset;
  }
  else {
    H5::DataSpace dataspace(dims.size(), dims.data());
    H5::DataSet dataset(file_->createDataSet(name.c_str(), type, dataspace));

    dataset.write(data, type, dataspace, H5P_DEFAULT);
    return dataset;
  }
}

void HDF5Writer::addAttribute(const H5::DataSet& set, const std::string& name,
                              const std::vector<hsize_t>& size, H5::DataType type, const void* data) {
  if (set.attrExists(name)) {
    return;
  }

  H5::DataSpace space(size.size(), size.data());
  auto attribute = set.createAttribute(name, type, space);
  attribute.write(type, data);
}

void HDF5Writer::addAttribute(const H5::DataSet& set, const std::string& name,
                              const std::string& value) {
  H5::StrType str_type(H5::PredType::C_S1, value.size());
  addAttribute(set, name, std::vector<hsize_t>{1}, str_type, value.c_str());
}

bool HDF5Writer::exists(const std::string& name) const {
  auto code = H5Gget_objinfo(file_id_, name.c_str(), 0, NULL);
  return code == 0;
}

}  // namespace io
}  // namespace dca
