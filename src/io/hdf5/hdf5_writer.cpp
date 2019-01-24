// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements hdf5_writer.hpp.

#include "dca/io/hdf5/hdf5_writer.hpp"

#include <fstream>
#include <stdexcept>

namespace dca {
namespace io {
// dca::io::

HDF5Writer::~HDF5Writer() {
  while (my_group.size())
    close_group();
  if (my_file != NULL)
    close_file();
}

bool HDF5Writer::fexists(const char* filename) {
  std::ifstream ifile(filename);
  return bool(ifile);
}

H5::H5File& HDF5Writer::open_file(std::string file_name, bool overwrite) {
  if (my_file != NULL or my_group.size() != 0)
    throw std::logic_error(__FUNCTION__);

  if (overwrite) {
    file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }
  else {
    if (fexists(file_name.c_str()))
      file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    else
      file_id = H5Fcreate(file_name.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  }

  if(file_id < 0)
    throw std::runtime_error("Cannot open file : " + file_name);

  my_file = new H5::H5File(file_name.c_str(), H5F_ACC_RDWR);

  return (*my_file);
}

void HDF5Writer::close_file() {
  delete my_file;
  my_file = NULL;

  H5Fclose(file_id);
}

void HDF5Writer::open_group(std::string name) {
  my_paths.push_back(name);
  my_group.push_back(NULL);

  my_group.back() = new H5::Group(my_file->createGroup(get_path().c_str()));
}

void HDF5Writer::close_group() {
  delete my_group.back();

  my_group.pop_back();
  my_paths.pop_back();
}

std::string HDF5Writer::get_path() {
  std::string path = "/";

  for (size_t i = 0; i < my_paths.size(); i++) {
    path = path + my_paths[i];

    if (i < my_paths.size() - 1)
      path = path + "/";
  }

  // cout << path << endl;

  return path;
}

void HDF5Writer::execute(const std::string& name,
                         const std::string& value)  //, H5File& file, std::string path)
{
  if (value.size() > 0) {
    H5::H5File& file = (*my_file);
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
    H5::H5File& file = (*my_file);

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

}  // io
}  // dca
