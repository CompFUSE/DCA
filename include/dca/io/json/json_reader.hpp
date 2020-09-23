// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// JSON reader.

#ifndef DCA_IO_JSON_JSON_READER_HPP
#define DCA_IO_JSON_JSON_READER_HPP

#include <iostream>
#include <stack>
#include <string>

#include "dca/io/json/details/json_group.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"

namespace dca::io {

class JSONReader {
public:
  JSONReader(bool verbose = true);
  ~JSONReader() = default;

  void open_file(const std::string& filename);
  void close_file() noexcept;

  void open_group(const std::string& name) noexcept;
  void close_group();

  constexpr static bool is_reader = true;
  constexpr static bool is_writer = false;

  template <class T>
  bool execute(const std::string& name, T& obj) noexcept;

  template <class Scalar, class Domain>
  bool execute(const std::string& name, func::function<Scalar, Domain>& f) noexcept;

  template <class Scalar>
  bool execute(const std::string& name, linalg::Matrix<Scalar, dca::linalg::CPU>& m);

  template <class T>
  bool execute(T& f) {
    return execute(f.get_name(), f);
  }

private:
  bool verbose_;
  std::stack<details::JSONGroup*> open_groups_;
  details::JSONGroup root_;
};

template <class T>
bool JSONReader::execute(const std::string& name, T& obj) noexcept {
  // TODO: perform check when opening group
  if (!open_groups_.top()) {
    return false;
  }

  return open_groups_.top()->readEntry(name, obj);
}

template <class Scalar, class Domain>
bool JSONReader::execute(const std::string& name, func::function<Scalar, Domain>& f) noexcept {
  if (verbose_)
    std::cout << "\t starts writing function : " << f.get_name() << "\n";

  open_group(name);

  bool result = true;

  std::vector<Scalar> data;
  std::vector<std::size_t> dmn_sizes;
  std::string f_name;

  result &= execute("name", f_name);
  result &= execute("domain-sizes", dmn_sizes);
  result &= execute("data", data);

  close_group();

  if (!result || dmn_sizes != f.getDomainSizes()) {
    if (verbose_) {
      std::cerr << "Domain sizes mismatch" << std::endl;
    }
    return false;
  }

  std::copy(data.begin(), data.end(), f.values());
  f.set_name(f_name);

  return true;
}

template <class Scalar>
bool JSONReader::execute(const std::string& name, linalg::Matrix<Scalar, dca::linalg::CPU>& m) {
  open_group(name);

  bool result = true;

  std::vector<std::vector<Scalar>> data;
  std::pair<int, int> size;
  std::string m_name;

  result &= execute("name", m_name);
  result &= execute("size", size);
  result &= execute("data", data);

  close_group();

  if (!result) {
    return false;
  }

  m.resize(size);

  for (int i = 0; i < size.first; ++i)
    for (int j = 0; j < size.second; ++j) {
      m(i, j) = data.at(i).at(j);
    }
  m.set_name(m_name);

  return true;
}

}  // namespace dca::io

#endif  // #ifndef DCA_IO_JSON_JSON_READER_HPP
