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

#include <fstream>
#include <stdexcept>

namespace dca {
namespace io {
// dca::io::

HDF5Writer::~HDF5Writer() {
  if (file_)
    close_file();
}

bool HDF5Writer::fexists(const char* filename) {
  std::ifstream ifile(filename);
  return bool(ifile);
}

void HDF5Writer::open_file(std::string file_name, bool overwrite) {
  if (file_)
    throw std::logic_error(__FUNCTION__);

  if (overwrite) {
    file_ = std::make_unique<H5::H5File>(file_name.c_str(), H5F_ACC_TRUNC);
  }
  else {
    if (fexists(file_name.c_str()))
      file_ = std::make_unique<H5::H5File>(file_name.c_str(), H5F_ACC_RDWR);
    else
      file_ = std::make_unique<H5::H5File>(file_name.c_str(), H5F_ACC_EXCL);
  }

  file_id_ = file_->getId();
}

void HDF5Writer::close_file() {
  file_->close();
  file_.release();
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
  std::string full_name = get_path() + '/' + name;

  write(full_name, std::vector<hsize_t>{value.size()}, HDF5_TYPE<char>::get_PredType(), value.data());
}

void HDF5Writer::execute(const std::string& name,
                         const std::vector<std::string>& value)  //, H5File& file, std::string path)
{
  if (value.size() > 0) {
    open_group(name);
    execute("size", value.size());

    open_group("data");

    const auto path = get_path();

    for (int i = 0; i < value.size(); ++i) {
      execute(get_path() + "/" + std::to_string(i), value[i]);
    }

    close_group();
    close_group();
  }
}

void HDF5Writer::write(const std::string& name, const std::vector<hsize_t>& dims, H5::PredType type,
                       const void* data) {
  if (exists(name)) {
    H5::DataSet dataset = file_->openDataSet(name.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

#ifndef NDEBUG  // Check for a size match.
    std::vector<hsize_t> size_check(dims.size());
    dataspace.getSimpleExtentDims(size_check.data(), nullptr);
    assert(dims == size_check);
#endif

    dataset.write(data, type, dataspace, H5P_DEFAULT);
  }
  else {
    H5::DataSpace dataspace(dims.size(), dims.data());
    H5::DataSet dataset(file_->createDataSet(name.c_str(), type, dataspace));

    dataset.write(data, type, dataspace, H5P_DEFAULT);
  }
}

bool HDF5Writer::exists(const std::string& name) {
  auto code = H5Gget_objinfo(file_id_, name.c_str(), 0, NULL);
  return code == 0;
}

}  // namespace io
}  // namespace dca
