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

bool fileExists(const std::string& filename);

class HDF5Writer {
public:
  typedef H5::H5File file_type;

public:
  // In: verbose. If true, the writer outputs a short log whenever it is executed.
  HDF5Writer(bool verbose = true) : verbose_(verbose) {
    H5::Exception::dontPrint();
  }

  ~HDF5Writer();

  constexpr bool is_reader() {
    return false;
  }
  constexpr bool is_writer() {
    return true;
  }

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

  template <typename Scalar>
  void execute(const std::string& name, const std::vector<std::complex<Scalar>>& value);

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
  void execute(const func::function<std::complex<Scalar>, domain_type>& f);

  template <typename Scalar, typename domain_type>
  void execute(const std::string& name, const func::function<Scalar, domain_type>& f);

  template <typename Scalar, typename domain_type>
  void execute(const std::string& name, const func::function<std::complex<Scalar>, domain_type>& f);

  template <typename Scalar>
  void execute(const std::string& name, const dca::linalg::Vector<Scalar, dca::linalg::CPU>& A);

  template <typename Scalar>
  void execute(const std::string& name,
               const dca::linalg::Vector<std::complex<Scalar>, dca::linalg::CPU>& A);

  template <typename Scalar>
  void execute(const std::string& name, const dca::linalg::Matrix<Scalar, dca::linalg::CPU>& A);

  template <typename Scalar>
  void execute(const std::string& name,
               const dca::linalg::Matrix<std::complex<Scalar>, dca::linalg::CPU>& A);

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

  operator bool() const {
    return static_cast<bool>(file_);
  }

  void lock() {
    mutex_.lock();
  }

  void unlock() {
    mutex_.unlock();
  }

  void set_verbose(bool verbose) {
    verbose_ = verbose;
  }

private:
  bool exists(const std::string& name) const;

  void write(const std::string& name, const std::vector<hsize_t>& size, H5::PredType type,
             const void* data);

  std::unique_ptr<H5::H5File> file_;

  hid_t file_id_;

  std::vector<std::string> my_paths_;

  bool verbose_;

  std::mutex mutex_;

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
void HDF5Writer::execute(const std::string& name, const std::vector<std::complex<Scalar>>& value) {
  if (value.size() > 0) {
    std::string full_name = get_path() + "/" + name;
    std::vector<hsize_t> dims{value.size(), 2};
    write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), value.data());
  }
}

template <typename Scalar>
void HDF5Writer::execute(const std::string& name, const std::vector<std::vector<Scalar>>& value) {
  if (value.size() > 0) {
    std::string full_name = get_path() + "/" + name;

    bool all_the_same_size = true;
    const std::size_t cols = value[0].size();
    for (auto& v : value) {
      if (v.size() != cols) {
        all_the_same_size = false;
        break;
      }
    }

    if (all_the_same_size) {
      std::vector<hsize_t> dims{value.size(), cols};
      std::vector<Scalar> linearized(dims[0] * dims[1]);

      for (std::size_t i = 0, linindex = 0; i < value.size(); ++i)
        for (std::size_t j = 0; j < cols; ++j)
          linearized[linindex++] = value[i][j];

      write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), linearized.data());
    }
    else {
      open_group(full_name);

      execute("size", value.size());

      open_group("data");

      std::vector<hsize_t> dims(1);
      for (std::size_t i = 0; i < value.size(); ++i) {
        const std::string new_name = full_name + "/data/row_" + std::to_string(i);
        dims[0] = value[i].size();
        write(new_name, dims, HDF5_TYPE<Scalar>::get_PredType(), value[i].data());
      }

      close_group();
      close_group();
    }
  }
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
void HDF5Writer::execute(const func::function<std::complex<Scalar>, domain_type>& f) {
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

  open_group(name);

  std::string new_path = get_path();

  execute("name", f.get_name());

  std::vector<hsize_t> dims;

  for (int l = 0; l < f.signature(); ++l)
    dims.push_back(f[l]);

  execute("domain-sizes", dims);

  // be carefull --> HDF5 is by default row-major, while the function-class is column-major !
  std::reverse(dims.begin(), dims.end());

  std::string full_name = new_path + "/data";

  write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), f.values());

  close_group();
}

template <typename Scalar, typename domain_type>
void HDF5Writer::execute(const std::string& name,
                         const func::function<std::complex<Scalar>, domain_type>& f) {
  if (f.size() == 0)
    return;

  open_group(name);

  std::string new_path = get_path();

  execute("name", f.get_name());

  std::vector<hsize_t> dims;

  for (int l = 0; l < f.signature(); ++l)
    dims.push_back(f[l]);

  execute("domain-sizes", dims);

  // be carefull --> HDF5 is by default row-major, while the function-class is column-major !
  std::reverse(dims.begin(), dims.end());

  dims.push_back(2);
  std::string full_name = get_path() + "/data";
  write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), f.values());

  close_group();
}

template <typename Scalar>
void HDF5Writer::execute(const std::string& name,
                         const dca::linalg::Vector<Scalar, dca::linalg::CPU>& V) {
  std::string full_name = get_path() + "/" + name;
  write(full_name, std::vector<hsize_t>{V.size()}, HDF5_TYPE<Scalar>::get_PredType(), V.ptr());
}

template <typename Scalar>
void HDF5Writer::execute(const std::string& name,
                         const dca::linalg::Vector<std::complex<Scalar>, dca::linalg::CPU>& V) {
  std::string full_name = get_path() + "/" + name;
  write(full_name, std::vector<hsize_t>{V.size(), 2}, HDF5_TYPE<Scalar>::get_PredType(), V.ptr());
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
  write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), linearized.data());
}

template <typename Scalar>
void HDF5Writer::execute(const std::string& name,
                         const dca::linalg::Matrix<std::complex<Scalar>, dca::linalg::CPU>& A) {
  std::vector<hsize_t> dims{hsize_t(A.nrRows()), hsize_t(A.nrCols()), hsize_t(2)};
  std::vector<std::complex<Scalar>> linearized(A.nrRows() * A.nrCols());

  int linindex = 0;
  // Note: Matrices are row major, while HDF5 is column major
  for (int i = 0; i < A.nrRows(); ++i)
    for (int j = 0; j < A.nrCols(); ++j)
      linearized[linindex++] = A(i, j);

  std::string full_name = get_path() + "/" + name;
  write(full_name, dims, HDF5_TYPE<Scalar>::get_PredType(), linearized.data());
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
