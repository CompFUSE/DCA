// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Norbert Podhorszki (pnorbert@ornl.gov)
//
// ADIOS2 reader.

#ifndef DCA_IO_ADIOS2_ADIOS2_READER_HPP
#define DCA_IO_ADIOS2_ADIOS2_READER_HPP

#include <complex>
#include <string>
#include <vector>

#include "adios2.h"

#include "dca/io/buffer.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace io {
// dca::io::

class ADIOS2Reader {
public:
  typedef adios2::ADIOS file_type;

public:
  // In: verbose. If true, the reader outputs a short log whenever it is executed.
  ADIOS2Reader(const std::string& config = "", bool verbose = true);
  ~ADIOS2Reader();

  constexpr bool is_reader() {
    return true;
  }
  constexpr bool is_writer() {
    return false;
  }

  void open_file(const std::string& file_name);
  void close_file();

  void open_group(const std::string& name);
  void close_group();

  std::string get_path(const std::string& name = "");

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

  template <typename Scalar>
  std::vector<size_t> getSize(const std::string& name);

  adios2::ADIOS adios_;
  const bool verbose_;

  adios2::IO io_;
  std::string io_name_;
  std::string file_name_;
  adios2::Engine file_;

  std::vector<std::string> my_paths_;
};

template <typename arbitrary_struct_t>
void ADIOS2Reader::from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name) {
  ADIOS2Reader reader_obj;
  reader_obj.open_file(file_name);
  arbitrary_struct.read_write(reader_obj);
  reader_obj.close_file();
}

template <typename Scalar>
bool ADIOS2Reader::execute(const std::string& name, Scalar& value) {
  std::string full_name = get_path(name);

  if (!exists(full_name)) {
    return false;
  }

  file_.Get<Scalar>(full_name, &value, adios2::Mode::Sync);
  return true;
}

template <typename Scalar>
bool ADIOS2Reader::execute(const std::string& name, std::vector<Scalar>& value) {
  std::string full_name = get_path(name);

  if (!exists(full_name)) {
    return false;
  }

  /* ADIOS will resize vector to match size of incoming data */
  file_.Get<Scalar>(full_name, value, adios2::Mode::Sync);
  return true;
}

template <typename Scalar>
bool ADIOS2Reader::execute(const std::string& name, std::vector<std::vector<Scalar>>& value) {
  throw(std::logic_error("ADIOS does not support vector of vectors. name = " + name + " (" +
                         __FUNCTION__ + ")"));
}

template <typename Scalar, std::size_t n>
bool ADIOS2Reader::execute(const std::string& name, std::vector<std::array<Scalar, n>>& value) {
  throw(std::logic_error("ADIOS does not support variable length arrays. name = " + name + " (" +
                         __FUNCTION__ + ")"));
}

template <typename Scalartype, typename domain_type>
bool ADIOS2Reader::execute(func::function<Scalartype, domain_type>& f) {
  return execute(f.get_name(), f);
}

template <typename Scalartype, typename domain_type>
bool ADIOS2Reader::execute(const std::string& name, func::function<Scalartype, domain_type>& f) {
  std::string full_name = get_path(name);

  if (!exists(full_name)) {
    std::cout << "\n\n\t the function (" + name + ") does not exist in path : " + get_path() +
                     "\n\n";
    return false;
  }

  std::cout << "\n\tstart ADIOS reading function : " << name << std::endl;

  const std::string sizeAttrName = full_name + "/domain-sizes";
  auto sizeAttr = io_.InquireAttribute<size_t>(sizeAttrName);

  if (sizeAttr) {
    try {
      // Read sizes.

      const std::vector<size_t>& dims = sizeAttr.Data();

      // Check sizes.
      if (sizeAttr.Data().size() != f.signature())
        throw(std::length_error("The number of domains is different"));
      for (int i = 0; i < f.signature(); ++i) {
        if (dims[i] != f[i])
          throw(std::length_error("The size of domain " + std::to_string(i) + " is different"));
      }
    }
    catch (std::length_error& err) {
      std::cerr << "Could not perform a size check on the function  " << name
                << " error: " << err.what() << std::endl;
    }

    file_.Get(full_name, f.values(), adios2::Mode::Sync);
  }
  else {
    std::cerr << "Could not perform a size check on the function  " << name << std::endl;
  }
  return true;
}

template <typename Scalar>
bool ADIOS2Reader::execute(const std::string& name, dca::linalg::Vector<Scalar, dca::linalg::CPU>& V) {
  std::string full_name = get_path(name);
  if (!exists(full_name)) {
    return false;
  }

  auto dims = getSize<Scalar>(full_name);
  assert(dims.size() == 1);
  V.resize(dims.at(0));

  file_.Get(full_name, V.ptr(), adios2::Mode::Sync);
  return true;
}

template <typename Scalar>
bool ADIOS2Reader::execute(const std::string& name, dca::linalg::Matrix<Scalar, dca::linalg::CPU>& A) {
  std::string full_name = get_path(name);
  if (!exists(full_name)) {
    return false;
  }

  auto dims = getSize<Scalar>(full_name);
  assert(dims.size() == 2);

  std::vector<Scalar> linearized(dims[0] * dims[1]);
  file_.Get(full_name, linearized.data(), adios2::Mode::Sync);

  // ADIOS2 is column major, while Matrix is row major.
  A.resizeNoCopy(std::make_pair(dims[0], dims[1]));
  for (int i = 0, linindex = 0; i < A.nrRows(); ++i) {
    for (int j = 0; j < A.nrCols(); ++j)
      A(i, j) = linearized[linindex++];
  }

  A.set_name(name);

  return true;
}

template <typename Scalar>
bool ADIOS2Reader::execute(dca::linalg::Matrix<Scalar, dca::linalg::CPU>& A) {
  return execute(A.get_name(), A);
}

template <typename Scalar>
std::vector<size_t> ADIOS2Reader::getSize(const std::string& name) {
  adios2::Variable<Scalar> var = io_.InquireVariable<Scalar>(name);
  if (var) {
    return var.Shape();
  }
  else {
    return std::vector<size_t>();
  }
}

}  // namespace io
}  // namespace dca

#endif  // DCA_IO_ADIOS2_ADIOS2_READER_HPP
