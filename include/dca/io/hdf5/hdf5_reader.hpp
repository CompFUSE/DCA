// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// HDF5 reader.

#ifndef DCA_IO_HDF5_HDF5_READER_HPP
#define DCA_IO_HDF5_HDF5_READER_HPP

#include <complex>
#include <string>
#include <vector>

#include "H5Cpp.h"

#include "dca/io/buffer.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_types.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace io {
// dca::io::

class HDF5Reader {
public:
  typedef H5::H5File file_type;

  HDF5Reader() : my_file(NULL), my_paths(0) {}
  ~HDF5Reader();

  bool is_reader() {
    return true;
  }
  bool is_writer() {
    return false;
  }

  void open_file(std::string file_name);
  void close_file();

  void open_group(std::string name) {
    my_paths.push_back(name);
  }
  void close_group() {
    my_paths.pop_back();
  }

  std::string get_path();

  template <typename arbitrary_struct_t>
  static void from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name);

  template <typename scalartype>
  void execute(std::string name, scalartype& value);

  template <typename scalar_type>
  void execute(std::string name, std::vector<scalar_type>& value);

  template <typename scalar_type>
  void execute(std::string name, std::vector<std::complex<scalar_type>>& value);

  void execute(std::string name, std::string& value);

  void execute(std::string name, std::vector<std::string>& value);

  // TODO: Remove? (only thing that depends on domains.hpp)
  template <typename domain_type>
  void execute(std::string /*name*/, func::dmn_0<domain_type>& /*dmn*/) {}

  template <typename scalartype, typename domain_type>
  void execute(func::function<scalartype, domain_type>& f);

  template <typename scalartype, typename domain_type>
  void execute(std::string name, func::function<scalartype, domain_type>& f);

  template <typename scalar_type>
  void execute(std::string name, dca::linalg::Vector<scalar_type, dca::linalg::CPU>& A);

  template <typename scalar_type>
  void execute(std::string name, dca::linalg::Vector<std::complex<scalar_type>, dca::linalg::CPU>& A);

  template <typename scalar_type>
  void execute(std::string name, dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A);

  void execute(std::string name, io::Buffer& buff) {
    return execute(name, static_cast<io::Buffer::Container&>(buff));
  }

private:
  bool fexists(const char* filename);

  H5::H5File* my_file;
  std::vector<std::string> my_paths;
};

template <typename arbitrary_struct_t>
void HDF5Reader::from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name) {
  HDF5Reader reader_obj;
  reader_obj.open_file(file_name);
  arbitrary_struct.read_write(reader_obj);
  reader_obj.close_file();
}

template <typename scalar_type>
void HDF5Reader::execute(std::string name, scalar_type& value) {
  std::string full_name = get_path() + "/" + name;

  try {
    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &value);
  }
  catch (...) {
    std::cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
  }
}

template <typename scalar_type>
void HDF5Reader::execute(std::string name, std::vector<scalar_type>& value) {
  std::string full_name = get_path() + "/" + name;

  try {
    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    value.resize(dataset.getInMemDataSize() / sizeof(scalar_type));

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &value[0]);
  }
  catch (...) {
    std::cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    // throw std::logic_error(__FUNCTION__);
  }
}

template <typename scalar_type>
void HDF5Reader::execute(std::string name, std::vector<std::complex<scalar_type>>& value) {
  std::string full_name = get_path() + "/" + name;

  try {
    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    value.resize(dataset.getInMemDataSize() / sizeof(std::complex<scalar_type>));

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &value[0]);
  }
  catch (...) {
    std::cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    // throw std::logic_error(__FUNCTION__);
  }
}

template <typename scalartype, typename domain_type>
void HDF5Reader::execute(func::function<scalartype, domain_type>& f) {
  execute(f.get_name(), f);
}

template <typename scalartype, typename domain_type>
void HDF5Reader::execute(std::string name, func::function<scalartype, domain_type>& f) {
  std::cout << "\n\tstart reading function : " << name;

  try {
    open_group(name);

    std::string full_name = get_path() + "/data";

    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalartype>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &f(0));

    close_group();
  }
  catch (...) {
    std::cout << "\n\n\t the function (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    // throw std::logic_error(__FUNCTION__);
  }
}

template <typename scalar_type>
void HDF5Reader::execute(std::string name, dca::linalg::Vector<scalar_type, dca::linalg::CPU>& V) {
  try {
    open_group(name);

    std::string full_name = get_path() + "/data";

    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    V.resize(dataset.getInMemDataSize() / sizeof(scalar_type));

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &V[0]);

    close_group();
  }
  catch (...) {
    std::cout << "\n\n\t the vector (" + name + ") does not exist in path : " + get_path() + "\n\n";
    // throw std::logic_error(__FUNCTION__);
  }
}

template <typename scalar_type>
void HDF5Reader::execute(std::string name,
                         dca::linalg::Vector<std::complex<scalar_type>, dca::linalg::CPU>& V) {
  try {
    open_group(name);

    std::string full_name = get_path() + "/data";

    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    V.resize(dataset.getInMemDataSize() / sizeof(std::complex<scalar_type>));

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &V[0]);

    close_group();
  }
  catch (...) {
    std::cout << "\n\n\t the vector (" + name + ") does not exist in path : " + get_path() + "\n\n";
  }
}

template <typename scalar_type>
void HDF5Reader::execute(std::string name, dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A) {
  try {
    open_group(name);

    std::string full_name = get_path() + "/data";

    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    H5::DataSpace dataspace = dataset.getSpace();

    // These 2 lines fix the bug of reading into a matrix which has been resized to a smaller size
    // hsize_t global_size[2] = {A.nrCols(), A.nrRows()}; // HDF5 use row
    // major data distribution
    // hsize_t global_size[2] = {A.capacity().second, A.capacity().first}; // HDF5 use
    // row major data distribution
    // dataspace.setExtentSimple(2, &global_size[0], NULL);

    H5Dread(dataset.getId(), HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &A(0, 0));

    close_group();
  }
  catch (...) {
    std::cout << "\n\n\t the function (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
  }
}

}  // io
}  // dca

#endif  // DCA_IO_HDF5_HDF5_READER_HPP
