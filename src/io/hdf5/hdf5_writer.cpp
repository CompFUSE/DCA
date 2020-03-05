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

H5::H5File& HDF5Writer::open_file(std::string file_name, bool overwrite) {
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

  file_id = file_->getId();

  if (file_id < 0)
    throw std::runtime_error("Cannot open file : " + file_name);

  return (*file_);
}

void HDF5Writer::close_file() {
  groups_.clear();
  file_.release();
}

void HDF5Writer::open_group(std::string name) {
  my_paths.push_back(name);
  const std::string path = get_path();

  if (!pathExists(path)) {
    groups_[path] = std::make_unique<H5::Group>(file_->createGroup(path.c_str()));
  }
}

void HDF5Writer::close_group() {
  my_paths.pop_back();
}

std::string HDF5Writer::get_path() {
  std::string path = "/";

  for (size_t i = 0; i < my_paths.size(); i++) {
    path += my_paths[i];

    if (i < my_paths.size() - 1)
      path += "/";
  }

  return path;
}

void HDF5Writer::execute(const std::string& name,
                         const std::string& value)  //, H5File& file, std::string path)
{
  if (value.size() > 0) {
    H5::H5File& file = (*file_);
    std::string path = get_path();

    hsize_t dims[1];

    H5::DataSet* dataset = NULL;
    H5::DataSpace* dataspace = NULL;

    {
      dims[0] = value.size();
      dataspace = new H5::DataSpace(1, dims);

      std::string full_name = path + "/" + name;
      dataset = new H5::DataSet(
          file.createDataSet(full_name.c_str(), HDF5_TYPE<char>::get_PredType(), *dataspace));

      H5Dwrite(dataset->getId(), HDF5_TYPE<char>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT,
               &value[0]);
    }

    delete dataset;
    delete dataspace;
  }
}

void HDF5Writer::execute(const std::string& name,
                         const std::vector<std::string>& value)  //, H5File& file, std::string path)
{
  if (value.size() > 0) {
    H5::H5File& file = (*file_);

    open_group(name);

    execute("size", value.size());  //, file, new_path);

    open_group("data");

    hsize_t dims[1];

    H5::DataSet* dataset = NULL;
    H5::DataSpace* dataspace = NULL;

    for (size_t l = 0; l < value.size(); l++) {
      dims[0] = value[l].size();
      dataspace = new H5::DataSpace(1, dims);

      std::stringstream ss;
      ss << get_path() << "/" << l;

      dataset = new H5::DataSet(
          file.createDataSet(ss.str().c_str(), HDF5_TYPE<char>::get_PredType(), *dataspace));

      H5Dwrite(dataset->getId(), HDF5_TYPE<char>::get(), dataspace->getId(), H5S_ALL, H5P_DEFAULT,
               &(value[l][0]));
    }

    close_group();

    delete dataset;
    delete dataspace;

    close_group();
  }
}

bool HDF5Writer::pathExists(const std::string& path) const {
  return groups_.count(path);
}

}  // namespace io
}  // namespace dca
