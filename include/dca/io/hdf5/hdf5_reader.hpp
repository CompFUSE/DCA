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

  // In: verbose. If true, the reader outputs a short log whenever it is executed.
  HDF5Reader(bool verbose = true) : my_file(NULL), my_paths(0), verbose_(verbose) {}

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

  // `execute` returns true if the object is read correctly.

  template <typename scalartype>
  bool execute(std::string name, scalartype& value);

  template <typename scalar_type>
  bool execute(std::string name, std::vector<scalar_type>& value);

  template <typename scalar_type>
  bool execute(std::string name, std::vector<std::complex<scalar_type>>& value);

  bool execute(std::string name, std::string& value);

  bool execute(std::string name, std::vector<std::string>& value);

  // TODO: Remove? (only thing that depends on domains.hpp)
  template <typename domain_type>
  bool execute(std::string /*name*/, func::dmn_0<domain_type>& /*dmn*/) {
    return false;
  }

  template <typename scalartype, typename domain_type>
  bool execute(func::function<scalartype, domain_type>& f);

  template <typename scalartype, typename domain_type>
  bool execute(std::string name, func::function<scalartype, domain_type>& f);

  template <typename scalar_type>
  bool execute(std::string name, dca::linalg::Vector<scalar_type, dca::linalg::CPU>& A);

  template <typename scalar_type>
  bool execute(std::string name, dca::linalg::Vector<std::complex<scalar_type>, dca::linalg::CPU>& A);

  template <typename scalar_type>
  bool execute(std::string name, dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A);

  template <typename scalar_type>
  bool execute(std::string name, dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU>& A);

  template <typename scalar_type>
  bool execute(dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A);

  bool execute(std::string name, io::Buffer& buff) {
    return execute(name, static_cast<io::Buffer::Container&>(buff));
  }

private:
  bool fexists(const char* filename);

  H5::H5File* my_file;
  std::vector<std::string> my_paths;

  bool verbose_;
};

template <typename arbitrary_struct_t>
void HDF5Reader::from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name) {
  HDF5Reader reader_obj;
  reader_obj.open_file(file_name);
  arbitrary_struct.read_write(reader_obj);
  reader_obj.close_file();
}

template <typename scalar_type>
bool HDF5Reader::execute(std::string name, scalar_type& value) {
  std::string full_name = get_path() + "/" + name;

  try {
    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &value);

    return true;
  }
  catch (...) {
    std::cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    return false;
  }
}

template <typename scalar_type>
bool HDF5Reader::execute(std::string name, std::vector<scalar_type>& value) {
  std::string full_name = get_path() + "/" + name;

  try {
    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    value.resize(dataset.getInMemDataSize() / sizeof(scalar_type));

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &value[0]);
    return true;
  }
  catch (...) {
    std::cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    return false;
  }
}

template <typename scalar_type>
bool HDF5Reader::execute(std::string name, std::vector<std::complex<scalar_type>>& value) {
  std::string full_name = get_path() + "/" + name;

  try {
    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    value.resize(dataset.getInMemDataSize() / sizeof(std::complex<scalar_type>));

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &value[0]);
    return true;
  }
  catch (...) {
    std::cout << "\n\n\t the variable (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    return false;
  }
}

template <typename scalartype, typename domain_type>
bool HDF5Reader::execute(func::function<scalartype, domain_type>& f) {
  return execute(f.get_name(), f);
}

template <typename scalartype, typename domain_type>
bool HDF5Reader::execute(std::string name, func::function<scalartype, domain_type>& f) {
  std::cout << "\n\tstart reading function : " << name;
  open_group(name);
  bool success = true;

  try {
    std::string full_name = get_path() + "/data";

    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalartype>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &f(0));
  }
  catch (const H5::FileIException& err) {
    std::cout << "\n\n\t the function (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    success = false;
  }

  close_group();
  return success;
}

template <typename scalar_type>
bool HDF5Reader::execute(std::string name, dca::linalg::Vector<scalar_type, dca::linalg::CPU>& V) {
  bool success = true;

  try {
    std::string full_name = get_path() + "/" + name;

    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    V.resize(dataset.getInMemDataSize() / sizeof(scalar_type));

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &V[0]);

    V.set_name(name);
  }
  catch (const H5::FileIException& err) {
    std::cout << "\n\n\t the vector (" + name + ") does not exist in path : " + get_path() + "\n\n";
    success = false;
  }

  return success;
}

template <typename scalar_type>
bool HDF5Reader::execute(std::string name,
                         dca::linalg::Vector<std::complex<scalar_type>, dca::linalg::CPU>& V) {
  bool success = true;

  try {
    std::string full_name = get_path() + "/" + name;

    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());

    V.resize(dataset.getInMemDataSize() / sizeof(std::complex<scalar_type>));

    H5::DataSpace dataspace = dataset.getSpace();

    H5Dread(dataset.getId(), HDF5_TYPE<scalar_type>::get(), dataspace.getId(), H5S_ALL, H5P_DEFAULT,
            &V[0]);

    V.set_name(name);
  }
  catch (const H5::FileIException& err) {
    std::cout << "\n\n\t the vector (" + name + ") does not exist in path : " + get_path() + "\n\n";
    success = false;
  }

  return success;
}

template <typename scalar_type>
bool HDF5Reader::execute(std::string name, dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A) {
  bool success = true;

  try {
    std::string full_name = get_path() + "/" + name;

    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    std::array<hsize_t, 2> dims;
    dataspace.getSimpleExtentDims(dims.data(), nullptr);

    std::vector<scalar_type> linearized(dims[0] * dims[1]);
    dataset.read(linearized.data(), HDF5_TYPE<scalar_type>::get_PredType(), dataspace);

    // HDF5 is column major, while Matrix is not major.
    A.resizeNoCopy(std::make_pair(dims[0], dims[1]));
    unsigned linindex = 0;
    for (int i = 0; i < A.nrRows(); ++i) {
      for (int j = 0; j < A.nrCols(); ++j)
        A(i, j) = linearized[linindex++];
    }

    A.set_name(name);
  }
  catch (const H5::FileIException& err) {
    std::cout << "\n\n\t the function (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    success = false;
  }

  return success;
}

template <typename scalar_type>
bool HDF5Reader::execute(std::string name,
                         dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU>& A) {
  bool success = true;

  try {
    std::string full_name = get_path() + "/" + name;

    H5::DataSet dataset = my_file->openDataSet(full_name.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    std::array<hsize_t, 3> dims;
    dataspace.getSimpleExtentDims(dims.data(), nullptr);

    A.resizeNoCopy(std::make_pair(dims[1], dims[2]));
    std::vector<std::complex<scalar_type>> linearized(A.nrRows() * A.nrCols());

    dataset.read(linearized.data(), HDF5_TYPE<scalar_type>::get_PredType(), dataspace);

    // HDF5 is column major, while Matrix is not major.
    unsigned linindex = 0;
    for (int i = 0; i < A.nrRows(); ++i)
      for (int j = 0; j < A.nrCols(); ++j) {
        A(i, j) = linearized[linindex++];
      }

    A.set_name(name);
  }
  catch (const H5::FileIException& err) {
    std::cout << "\n\n\t the function (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    success = false;
  }

  return success;
}

template <typename scalar_type>
bool HDF5Reader::execute(dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A) {
  return execute(A.get_name(), A);
}

}  // namespace io
}  // namespace dca

#endif  // DCA_IO_HDF5_HDF5_READER_HPP
