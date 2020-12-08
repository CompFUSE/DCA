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
  HDF5Reader(bool verbose = true) : verbose_(verbose) {}

  ~HDF5Reader();

  constexpr static bool is_reader = true;
  constexpr static bool is_writer = false;

  void open_file(std::string file_name);
  void close_file();

  void open_group(std::string name) {
    paths_.push_back(name);
  }
  void close_group() {
    paths_.pop_back();
  }

  std::string get_path();

  template <typename arbitrary_struct_t>
  static void from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name);

  // `execute` returns true if the object is read correctly.

  template <typename Scalartype>
  bool execute(const std::string& name, Scalartype& value);

  template <typename Scalar>
  bool execute(const std::string& name, std::vector<Scalar>& value);

  template <typename Scalar>
  bool execute(const std::string& name, std::vector<std::vector<Scalar>>& value);

  template <typename Scalar, std::size_t n>
  bool execute(const std::string& name, std::vector<std::array<Scalar, n>>& value);

  bool execute(const std::string& name, std::string& value);

  bool execute(const std::string& name, std::vector<std::string>& value);

  // TODO: Remove? (only thing that depends on domains.hpp)
  template <typename domain_type>
  bool execute(std::string /*name*/, func::dmn_0<domain_type>& /*dmn*/) {
    return false;
  }

  template <typename Scalartype, typename domain_type>
  bool execute(func::function<Scalartype, domain_type>& f);

  template <typename Scalartype, typename domain_type>
  bool execute(const std::string& name, func::function<Scalartype, domain_type>& f);

  template <typename Scalar>
  bool execute(const std::string& name, dca::linalg::Vector<Scalar, dca::linalg::CPU>& A);

  template <typename Scalar>
  bool execute(const std::string& name, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& A);

  template <typename Scalar>
  bool execute(dca::linalg::Matrix<Scalar, dca::linalg::CPU>& A);

  bool execute(const std::string& name, io::Buffer& buff) {
    return execute(name, static_cast<io::Buffer::Container&>(buff));
  }

private:
  bool exists(const std::string& name) const;

  void read(const std::string& name, H5::DataType type, void* data) const;
  std::vector<hsize_t> readSize(const std::string& name) const;

  std::unique_ptr<H5::H5File> file_;
  std::vector<std::string> paths_;

  bool verbose_;
};

template <typename arbitrary_struct_t>
void HDF5Reader::from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name) {
  HDF5Reader reader_obj;
  reader_obj.open_file(file_name);
  arbitrary_struct.read_write(reader_obj);
  reader_obj.close_file();
}

template <typename Scalar>
bool HDF5Reader::execute(const std::string& name, Scalar& value) {
  std::string full_name = get_path() + "/" + name;

  if (!exists(full_name)) {
    return false;
  }

  read(full_name, HDF5_TYPE<Scalar>::get_PredType(), &value);
  return true;
}

template <typename Scalar>
bool HDF5Reader::execute(const std::string& name, std::vector<Scalar>& value) {
  std::string full_name = get_path() + "/" + name;

  if (!exists(full_name)) {
    return false;
  }

  auto dims = readSize(full_name);
  assert(dims.size() == 1);
  value.resize(dims.at(0));

  read(full_name, HDF5_TYPE<Scalar>::get_PredType(), value.data());
  return true;
}

template <typename Scalar>
bool HDF5Reader::execute(const std::string& name, std::vector<std::vector<Scalar>>& value) {
  std::string full_name = get_path() + "/" + name;
  if (!exists(full_name)) {
    return false;
  }

  auto size = readSize(full_name)[0];
  const auto type = H5::VarLenType(HDF5_TYPE<Scalar>::get_PredType());

  std::vector<hvl_t> data(size);

  H5::DataSet dataset = file_->openDataSet(name.c_str());
  dataset.read(data.data(), type);

  value.resize(size);
  for (int i = 0; i < size; ++i) {
    value[i].resize(data[i].len);
    std::copy_n(static_cast<Scalar*>(data[i].p), data[i].len, value[i].data());
  }

  dataset.vlenReclaim(data.data(), type, dataset.getSpace());

  return true;
}

template <typename Scalar, std::size_t n>
bool HDF5Reader::execute(const std::string& name, std::vector<std::array<Scalar, n>>& value) {
  std::string full_name = get_path() + "/" + name;
  if (!exists(full_name)) {
    return false;
  }

  auto dims = readSize(full_name);
  assert(dims.size() == 2);
  if (dims.at(1) != n) {
    throw(std::length_error("Wrong array size"));
  }

  value.resize(dims[0]);
  read(full_name, HDF5_TYPE<Scalar>::get_PredType(), value.data());

  return true;
}

template <typename Scalartype, typename domain_type>
bool HDF5Reader::execute(func::function<Scalartype, domain_type>& f) {
  return execute(f.get_name(), f);
}

template <typename Scalartype, typename domain_type>
bool HDF5Reader::execute(const std::string& name, func::function<Scalartype, domain_type>& f) {
  std::string full_name = get_path() + "/" + name;

  if (!exists(full_name)) {
    std::cout << "\n\n\t the function (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    return false;
  }

  std::cout << "\n\tstart reading function : " << name;

  H5::DataSet dataset = file_->openDataSet(full_name.c_str());

  try {
    // Read sizes.
    std::vector<hsize_t> dims;
    auto domain_attribute = dataset.openAttribute("domain-sizes");
    hsize_t n_dims;
    domain_attribute.getSpace().getSimpleExtentDims(&n_dims);
    dims.resize(n_dims);
    domain_attribute.read(HDF5_TYPE<hsize_t>::get_PredType(), dims.data());

    // Check sizes.
    if (dims.size() != f.signature())
      throw(std::length_error("The number of domains is different"));
    for (int i = 0; i < f.signature(); ++i) {
      if (dims[i] != f[i])
        throw(std::length_error("The size of domain " + std::to_string(i) + " is different"));
    }
  }
  catch (H5::Exception& err) {
    std::cerr << "Could not perform a size check on the function  " << name << std::endl;
  }

  read(full_name, HDF5_TYPE<Scalartype>::get_PredType(), f.values());

  return true;
}

template <typename Scalar>
bool HDF5Reader::execute(const std::string& name, dca::linalg::Vector<Scalar, dca::linalg::CPU>& V) {
  std::string full_name = get_path() + "/" + name;
  if (!exists(full_name)) {
    return false;
  }

  auto dims = readSize(full_name);
  assert(dims.size() == 1);
  V.resize(dims.at(0));

  read(full_name, HDF5_TYPE<Scalar>::get_PredType(), V.ptr());

  return true;
}

template <typename Scalar>
bool HDF5Reader::execute(const std::string& name, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& A) {
  std::string full_name = get_path() + "/" + name;
  if (!exists(full_name)) {
    return false;
  }

  auto dims = readSize(full_name);
  assert(dims.size() == 2);

  std::vector<Scalar> linearized(dims[0] * dims[1]);
  read(full_name, HDF5_TYPE<Scalar>::get_PredType(), linearized.data());

  // HDF5 is column major, while Matrix is row major.
  A.resizeNoCopy(std::make_pair(dims[0], dims[1]));
  for (int i = 0, linindex = 0; i < A.nrRows(); ++i) {
    for (int j = 0; j < A.nrCols(); ++j)
      A(i, j) = linearized[linindex++];
  }

  A.set_name(name);

  return true;
}

template <typename Scalar>
bool HDF5Reader::execute(dca::linalg::Matrix<Scalar, dca::linalg::CPU>& A) {
  return execute(A.get_name(), A);
}

}  // namespace io
}  // namespace dca

#endif  // DCA_IO_HDF5_HDF5_READER_HPP
