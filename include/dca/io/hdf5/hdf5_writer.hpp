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
// HDF5 writer.

#ifndef DCA_IO_HDF5_HDF5_WRITER_HPP
#define DCA_IO_HDF5_HDF5_WRITER_HPP

#include <complex>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include "H5Cpp.h"

#include "dca/function/domains.hpp"
#include "dca/io/buffer.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_types.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace io {
// dca::io

class HDF5Writer {
public:
  typedef H5::H5File file_type;

public:
  // In: verbose. If true, the writer outputs a short log whenever it is executed.
  HDF5Writer(bool verbose = true) : verbose_(verbose) {
    H5::Exception::dontPrint();
  }

  ~HDF5Writer();

  constexpr static bool is_reader = false;
  constexpr static bool is_writer = true;

  void open_file(std::string file_name_ref, bool overwrite = true);
  void close_file();

  void open_group(std::string new_path);
  void close_group();

  std::string get_path();

  template <typename arbitrary_struct_t>
  static void to_file(const arbitrary_struct_t& arbitrary_struct, const std::string& file_name);

  template <typename Scalar>
  void execute(const std::string& name, Scalar value);

  template <typename Scalar>
  void execute(const std::string& name, const std::pair<Scalar, Scalar>& value);

  template <typename Scalar>
  void execute(const std::string& name, const std::vector<Scalar>& value);

  void execute(const std::string& name, const std::string& value);

  void execute(const std::string& name, const std::vector<std::string>& value);

  template <typename Scalar, std::size_t n>
  void execute(const std::string& name, const std::vector<std::array<Scalar, n>>& value);

  template <typename Scalar>
  void execute(const std::string& name, const std::vector<std::vector<Scalar>>& value);

  template <typename domain_type>
  void execute(const std::string& name, const func::dmn_0<domain_type>& dmn);

  template <typename Scalar, typename domain_type>
  void execute(const func::function<Scalar, domain_type>& f);

  template <typename Scalar, typename domain_type>
  void execute(const std::string& name, const func::function<Scalar, domain_type>& f);

  template <typename Scalar>
  void execute(const std::string& name, const dca::linalg::Vector<Scalar, dca::linalg::CPU>& A);

  template <typename Scalar>
  void execute(const std::string& name, const dca::linalg::Matrix<Scalar, dca::linalg::CPU>& A);

  template <typename Scalar>
  void execute(const dca::linalg::Matrix<Scalar, dca::linalg::CPU>& A) {
    execute(A.get_name(), A);
  }

  template <class T>
  void execute(const std::string& name, const std::unique_ptr<T>& obj);

  template <class T>
  void execute(const std::unique_ptr<T>& obj);

  void execute(const std::string& name, const io::Buffer& buffer) {
    return execute(name, static_cast<io::Buffer::Container>(buffer));
  }

  operator bool() const noexcept {
    return static_cast<bool>(file_);
  }

  void set_verbose(bool verbose) {
    verbose_ = verbose;
  }

private:
  bool exists(const std::string& name) const;

  H5::DataSet write(const std::string& name, const std::vector<hsize_t>& size, H5::DataType type,
                    const void* data);
  void addAttribute(const H5::DataSet& set, const std::string& name,
                    const std::vector<hsize_t>& size, H5::DataType type, const void* data);
  void addAttribute(const H5::DataSet& set, const std::string& name, const std::string& value);

  std::unique_ptr<H5::H5File> file_;

  hid_t file_id_;

  std::vector<std::string> my_paths_;

  bool verbose_;

  std::vector<hsize_t> size_check_;
};

template <typename arbitrary_struct_t>
void HDF5Writer::to_file(const arbitrary_struct_t& arbitrary_struct, const std::string& file_name) {
  HDF5Writer wr_obj;
  wr_obj.open_file(file_name);
  arbitrary_struct.read_write(wr_obj);
  wr_obj.close_file();
}

template <typename Scalar>
void HDF5Writer::execute(const std::string& name, Scalar value) {
  const std::string full_name = get_path() + "/" + name;
  std::vector<hsize_t> dims{1};

  write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), &value);
}

template <typename Scalar>
void HDF5Writer::execute(const std::string& name, const std::pair<Scalar, Scalar>& value) {
  std::string full_name = get_path() + "/" + name;
  std::vector<hsize_t> dims{2};

  write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), &value.first);
}

template <typename Scalar>
void HDF5Writer::execute(const std::string& name,
                         const std::vector<Scalar>& value)  //, H5File& file, std::string path)
{
  if (value.size() > 0) {
    std::string full_name = get_path() + "/" + name;
    std::vector<hsize_t> dims{value.size()};
    write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), value.data());
  }
}

template <typename Scalar>
void HDF5Writer::execute(const std::string& name, const std::vector<std::vector<Scalar>>& value) {
  std::string full_name = get_path() + "/" + name;

  std::vector<hvl_t> data(value.size());
  for (int i = 0; i < value.size(); ++i) {
    data[i].p = const_cast<void*>(static_cast<const void*>((value[i].data())));
    data[i].len = value[i].size();
  }

  const auto type = H5::VarLenType(HDF5_TYPE<Scalar>::get_PredType());

  write(full_name, std::vector<hsize_t>{data.size()}, type, data.data());
}

template <typename Scalar, std::size_t n>
void HDF5Writer::execute(const std::string& name, const std::vector<std::array<Scalar, n>>& value) {
  if (value.size() == 0)
    return;

  std::vector<hsize_t> dims{value.size(), n};
  std::string full_name = get_path() + "/" + name;

  write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), value.data());
}

template <typename domain_type>
void HDF5Writer::execute(const std::string& name, const func::dmn_0<domain_type>& dmn) {
  open_group(name);

  execute("name", dmn.get_name());
  execute("elements", dmn.get_elements());

  close_group();
}

template <typename Scalar, typename domain_type>
void HDF5Writer::execute(const func::function<Scalar, domain_type>& f) {
  if (f.size() == 0)
    return;

  if (verbose_)
    std::cout << "\t starts writing function : " << f.get_name() << "\n";

  execute(f.get_name(), f);
}

template <typename Scalar, typename domain_type>
void HDF5Writer::execute(const std::string& name, const func::function<Scalar, domain_type>& f) {
  if (f.size() == 0)
    return;

  const std::string full_name = get_path() + "/" + name;

  std::vector<hsize_t> dims;
  for (int l = 0; l < f.signature(); ++l)
    dims.push_back(f[l]);

  // be carefull --> HDF5 is by default row-major, while the function-class is column-major !
  std::reverse(dims.begin(), dims.end());

  auto dataset = write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), f.values());

  addAttribute(dataset, "name", f.get_name());

  std::reverse(dims.begin(), dims.end());
  auto type = HDF5_TYPE<hsize_t>::get_PredType();
  addAttribute(dataset, "domain-sizes", std::vector<hsize_t>{dims.size()}, type, dims.data());
}

template <typename Scalar>
void HDF5Writer::execute(const std::string& name,
                         const dca::linalg::Vector<Scalar, dca::linalg::CPU>& V) {
  std::string full_name = get_path() + "/" + name;
  auto dataset =
      write(full_name, std::vector<hsize_t>{V.size()}, HDF5_TYPE<Scalar>::get_PredType(), V.ptr());

  addAttribute(dataset, "name", V.get_name());
}

template <typename Scalar>
void HDF5Writer::execute(const std::string& name,
                         const dca::linalg::Matrix<Scalar, dca::linalg::CPU>& A) {
  std::vector<hsize_t> dims{hsize_t(A.nrRows()), hsize_t(A.nrCols())};
  std::vector<Scalar> linearized(dims[0] * dims[1]);

  int linindex = 0;
  // Note: Matrices are row major, while HDF5 is column major
  for (int i = 0; i < A.nrRows(); ++i)
    for (int j = 0; j < A.nrCols(); ++j)
      linearized[linindex++] = A(i, j);

  std::string full_name = get_path() + "/" + name;
  auto dataset = write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), linearized.data());

  addAttribute(dataset, "name", A.get_name());
}

template <class T>
void HDF5Writer::execute(const std::string& name, const std::unique_ptr<T>& obj) {
  if (obj)
    execute(name, *obj);
}

template <class T>
void HDF5Writer::execute(const std::unique_ptr<T>& obj) {
  if (obj)
    execute(*obj);
}

}  // namespace io
}  // namespace dca

#endif  // DCA_IO_HDF5_HDF5_WRITER_HPP
