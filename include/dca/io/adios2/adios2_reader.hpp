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
#include <numeric>
#include <string>
#include <sstream>
#include <vector>

#include "adios2.h"

#include "dca/io/buffer.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

#include "dca/config/haves_defines.hpp"
#ifdef DCA_HAVE_MPI
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
#endif

namespace dca {
namespace io {
// dca::io::

class ADIOS2Reader {
public:
  typedef adios2::ADIOS file_type;

public:
  // In: verbose. If true, the reader outputs a short log whenever it is executed.
  ADIOS2Reader(const std::string& config = "", bool verbose = false);
#ifdef DCA_HAVE_MPI
  ADIOS2Reader(const dca::parallel::MPIConcurrency* concurrency, const std::string& config = "",
               bool verbose = false);
#endif

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

  template <typename Scalartype, typename domain_type, DistType DT>
  bool execute(func::function<Scalartype, domain_type, DT>& f);

  template <typename Scalartype, typename domain_type, DistType DT>
  bool execute(func::function<Scalartype, domain_type, DT>& f, uint64_t start, uint64_t end);

  template <typename Scalartype, typename domain_type, DistType DT>
  bool execute(func::function<Scalartype, domain_type, DT>& f, const std::vector<int>& start,
               const std::vector<int>& end);

  template <typename Scalartype, typename domain_type, DistType DT>
  bool execute(const std::string& name, func::function<Scalartype, domain_type, DT>& f);

  template <typename Scalartype, typename domain_type, DistType DT>
  bool execute(const std::string& name, func::function<Scalartype, domain_type, DT>& f, uint64_t start,
               uint64_t end);

  template <typename Scalartype, typename domain_type, DistType DT>
  bool execute(const std::string& name, func::function<Scalartype, domain_type, DT>& f,
               const std::vector<int>& start, const std::vector<int>& end);

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

  template <class T>
  std::string VectorToString(const std::vector<T>& v);

  adios2::ADIOS adios_;
  const bool verbose_;
#ifdef DCA_HAVE_MPI
  const dca::parallel::MPIConcurrency* concurrency_;
#endif

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
  std::string full_name = get_path(name);

  adios2::Variable<Scalar> var = io_.InquireVariable<Scalar>(full_name);
  if (!var) {
    return false;
  }

  auto aSizes = io_.InquireAttribute<size_t>("_vector_sizes", full_name);
  if (!aSizes) {
    throw(std::runtime_error("ADIOS2Reader: Variable " + full_name + " is not a vector of vectors\n"));
  }

  if (var.Shape().size() != 1) {
    throw(std::runtime_error(
        "ADIOS2Reader: Variable " + full_name +
        " is not a vector of vectors because it is not stored as a 1D array but as " +
        std::to_string(var.Shape().size()) + "\n"));
  }
  size_t nTotal = var.Shape()[0];
  std::vector<size_t> nSizes = aSizes.Data();
  size_t nVectors = nSizes.size();
  size_t sum = std::accumulate(nSizes.begin(), nSizes.end(), 0);

  if (sum != nTotal) {
    throw(std::runtime_error("ADIOS2Reader: Variable " + full_name +
                             " vector of vectors data corrupted since size of array is " +
                             std::to_string(nTotal) + " but sum of components is " +
                             std::to_string(sum) + "\n"));
  }

  value.resize(nVectors);

  if (verbose_) {
    std::cout << "\t ADIOS2Reader: Read Vector of " << nVectors << " vectors " << full_name
              << " sizes (";
  }

  size_t startpos = 0;
  for (int i = 0; i < nVectors; ++i) {
    size_t vsize = nSizes[i];
    if (verbose_) {
      std::cout << vsize;
      if (i < nVectors - 1) {
        std::cout << ", ";
      }
    }
    std::vector<Scalar>& vec = value[i];
    vec.resize(vsize);
    var.SetSelection({{startpos}, {vsize}});
    file_.Get(full_name, vec.data(), adios2::Mode::Sync);
    startpos += vsize;
  }
  if (verbose_) {
    std::cout << ")\n";
  }
  return true;
}

template <typename Scalar, std::size_t n>
bool ADIOS2Reader::execute(const std::string& name, std::vector<std::array<Scalar, n>>& value) {
  std::string full_name = get_path(name);
  if (!exists(full_name)) {
    return false;
  }
  auto dims = getSize<Scalar>(full_name);
  if (dims.size() == 0) {
    return false;
  }
  assert(dims.size() == 2);
  if (dims.at(1) != n) {
    throw(std::length_error("Wrong array size " + std::to_string(n) +
                            " when reading vector of array variable " + full_name +
                            " which has size {" + std::to_string(dims[0]) + "," +
                            std::to_string(dims[1]) + "}\n"));
  }
  value.resize(dims.at(0));

  if (verbose_) {
    std::cout << "\t ADIOS2Reader: Vector of Array " << full_name << " size {" << dims[0] << ", "
              << n << "} -> value {" << value.size() << ", " << n << "}\n";
  }

  const auto beginptr = &value[0][0];
  const auto endptr = &value[dims[0] - 1][n - 1] + 1;

  if (verbose_) {
    std::cout << "\t ADIOS2Reader: Check if input is contiguous. Target Array beginptr = " << beginptr
              << " endptr = " << endptr << " diff = " << endptr - beginptr
              << " expected = " << dims[0] * n << "\n";
  }

  assert(endptr - beginptr == dims[0] * n);

  Scalar* ptr = value[0].data();
  file_.Get(full_name, ptr, adios2::Mode::Sync);
  return true;
}

template <typename Scalartype, typename domain_type, DistType DT>
bool ADIOS2Reader::execute(func::function<Scalartype, domain_type, DT>& f) {
  return execute(f.get_name(), f);
}

template <typename Scalartype, typename domain_type, DistType DT>
bool ADIOS2Reader::execute(const std::string& name, func::function<Scalartype, domain_type, DT>& f) {
  std::string full_name = get_path(name);

  if (!exists(full_name)) {
    std::cout << "\n\t ADIOS2Reader:the function (" + name +
                     ") does not exist in path : " + get_path() + "\n\n";
    return false;
  }

  if (verbose_) {
    std::cout << "\t ADIOS2Reader: read function : " << name << std::endl;
  }

  const std::string sizeAttrName = full_name + "/domain-sizes";
  auto sizeAttr = io_.InquireAttribute<size_t>(sizeAttrName);

  if constexpr ( DT == dca::DistType::BLOCKED ) {
    execute(name, f, f.get_start_subindex(), f.get_end_subindex());
  }
  else {
  if (sizeAttr) {
    try {
      // Read sizes.

      const std::vector<size_t>& dims = sizeAttr.Data();

      // Check sizes.
      if (sizeAttr.Data().size() != f.signature())
        throw(std::length_error("The number of domains is different"));
      for (int i = 0; i < f.signature(); ++i) {
        if (dims[i] != f.get_domain().get_subdomain_size(i))
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
  }
  return true;
}

template <typename Scalartype, typename domain_type, DistType DT>
bool ADIOS2Reader::execute(func::function<Scalartype, domain_type, DT>& f, uint64_t start, uint64_t end) {
  return execute(f.get_name(), f, start, end);
}

template <typename Scalartype, typename domain_type, DistType DT>
bool ADIOS2Reader::execute(const std::string& name, func::function<Scalartype, domain_type, DT>& f,
                           uint64_t start, uint64_t end) {
  std::string full_name = get_path(name);
  adios2::Variable<Scalartype> var = io_.InquireVariable<Scalartype>(full_name);

  if (!var) {
    std::cout << "\n\t ADIOS2Reader:the function (" + name +
                     ") does not exist in path : " + get_path() + "\n\n";
    return false;
  }

  const int ndim = f.signature();

  if (var.Shape().size() != 1) {
    if (var.Shape().size() == ndim) {
      // get the N-dimensional decomposition
      std::vector<int> subind_start = f.linind_2_subind(start);
      std::vector<int> subind_end = f.linind_2_subind(end);
      return execute(name, f, subind_start, subind_end);
    }
    else {
      std::cerr << "ERROR ADIOS2Reader: Dimension check failed on function " << name
                << ". Reason: Function in file is not 1D but " << std::to_string(var.Shape().size())
                << " while the function in memory is " << std::to_string(ndim) << std::endl;
    }
  }

  if (verbose_) {
    std::cout << "\t ADIOS2Reader: Read function : " << name
              << " in linear distributed manner, rank = " << concurrency_->id()
              << " start = " << std::to_string(start) << " end = " << std::to_string(end) << "\n";
  }

  const std::string sizeAttrName = full_name + "/domain-sizes";
  auto sizeAttr = io_.InquireAttribute<size_t>(sizeAttrName);

  if (sizeAttr) {
    try {
      // Read sizes.
      const std::vector<size_t>& dims = sizeAttr.Data();

      // Check sizes.
      if (dims.size() != ndim)
        throw(std::length_error("The number of domains is different"));
      for (int i = 0; i < ndim; ++i) {
        if (dims[i] != f.get_domain().get_subdomain_size(i))
          throw(std::length_error("The size of domain " + std::to_string(i) + " is different"));
      }
    }
    catch (std::length_error& err) {
      std::cerr << "ERROR ADIOS2Reader: Size check failed on function " << name
                << ". Reason: " << err.what() << std::endl;
      return false;
    }

    // Make 1-D selection
    // ADIOS2 takes size_t vector of start and count
    // and in row-major (reverse) order
    std::vector<size_t> s = {start};
    std::vector<size_t> c = {end - start + 1};

    if (verbose_) {
      std::cout << "\t ADIOS2Reader: Read function : " << name
                << " in linear distributed manner, rank = " << concurrency_->id()
                << " shape = " << VectorToString(var.Shape()) << " start = " << VectorToString(s)
                << " count = " << VectorToString(c) << "\n";
    }

    var.SetSelection({s, c});
    /* TODO: must pass the correct data pointer here, not f.values() */
    file_.Get(full_name, f.values(), adios2::Mode::Sync);
  }
  else {
    std::cerr << "ERROR ADIOS2Reader: Size check failed on function " << name
              << ". Reason: missing attribute " << sizeAttrName << std::endl;
    return false;
  }
  return true;
}

template <typename Scalartype, typename domain_type, DistType DT>
bool ADIOS2Reader::execute(func::function<Scalartype, domain_type, DT>& f,
                           const std::vector<int>& start, const std::vector<int>& end) {
  return execute(f.get_name(), f, start, end);
}

template <typename Scalartype, typename domain_type, DistType DT>
bool ADIOS2Reader::execute(const std::string& name, func::function<Scalartype, domain_type, DT>& f,
                           const std::vector<int>& start, const std::vector<int>& end) {
  std::string full_name = get_path(name);
  adios2::Variable<Scalartype> var = io_.InquireVariable<Scalartype>(full_name);

  if (!var) {
    std::cout << "\n\t ADIOS2Reader:the function (" + name +
                     ") does not exist in path : " + get_path() + "\n\n";
    return false;
  }

  const int ndim = f.signature();

  if (start.size() != ndim || end.size() != ndim) {
    std::cerr << "ERROR ADIOS2Reader: Dimension mismatch on function " << name << ": function is "
              << std::to_string(ndim) << "-D, start is " << std::to_string(start.size())
              << "-D and end is " << std::to_string(end.size()) << std::endl;
    return false;
  }

  if (verbose_) {
    std::cout << "\t ADIOS2Reader: Read function : " << name
              << " in distributed manner, rank = " << concurrency_->id()
              << " start = " << VectorToString(start) << " end = " << VectorToString(end) << "\n";
  }

  const std::string sizeAttrName = full_name + "/domain-sizes";
  auto sizeAttr = io_.InquireAttribute<size_t>(sizeAttrName);

  if (sizeAttr) {
    try {
      // Read sizes.
      const std::vector<size_t>& dims = sizeAttr.Data();

      if (dims.size() != var.Shape().size()) {
        throw(std::length_error("The function " + name + " in the file was written as a " +
                                std::to_string(var.Shape().size()) +
                                "-D array. Cannot read it back with " +
                                std::to_string(start.size()) + "-D read selections"));
      }

      // Check sizes.
      if (dims.size() != f.signature())
        throw(std::length_error("The number of domains is different"));
      for (int i = 0; i < f.signature(); ++i) {
        if (dims[i] != f.get_domain().get_subdomain_size(i))
          throw(std::length_error("The size of domain " + std::to_string(i) + " is different"));
      }
    }
    catch (std::length_error& err) {
      std::cerr << "ERROR ADIOS2Reader: Size check failed on function " << name
                << ". Reason: " << err.what() << std::endl;
      return false;
    }

    // Make N-D selection
    // ADIOS2 takes size_t vector of start and count
    // and in row-major (reverse) order
    std::vector<size_t> s(ndim);
    std::vector<size_t> c(ndim);
    for (int i = 0; i < ndim; ++i) {
      s[ndim - i - 1] = static_cast<size_t>(start[i]);
      c[ndim - i - 1] = static_cast<size_t>(end[i] - start[i] + 1);
    }

    if (verbose_) {
      std::cout << "\t ADIOS2Reader: Read function : " << name
                << " in distributed manner, rank = " << concurrency_->id()
                << " shape = " << VectorToString(var.Shape()) << " start = " << VectorToString(s)
                << " count = " << VectorToString(c) << "\n";
    }

    var.SetSelection({s, c});
    /* TODO: must pass the correct data pointer here, not f.values() */
    file_.Get(full_name, f.values(), adios2::Mode::Sync);
  }
  else {
    std::cerr << "ERROR ADIOS2Reader: Size check failed on function " << name
              << ". Reason: missing attribute " << sizeAttrName << std::endl;
    return false;
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

template <class T>
std::string ADIOS2Reader::VectorToString(const std::vector<T>& v) {
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < v.size(); ++i) {
    if (i != 0)
      ss << ",";
    ss << v[i];
  }
  ss << "]";
  return ss.str();
}

}  // namespace io
}  // namespace dca

#endif  // DCA_IO_ADIOS2_ADIOS2_READER_HPP
