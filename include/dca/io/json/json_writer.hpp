// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// JSON writer.

#ifndef DCA_IO_JSON_JSON_WRITER_HPP
#define DCA_IO_JSON_JSON_WRITER_HPP

#include <fstream>
#include <string>
#include <stack>

#include "dca/io/json/details/json_group.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"

namespace dca::io {

class JSONWriter {
public:
  JSONWriter(bool verbose = true);
  ~JSONWriter();

  // Acquires lock on file. If "overwrite" is true and the file exists, the content is read into
  // memory and serves as the initial state of the writer.
  void open_file(const std::string& filename, bool overwrite = true);

  // Flushes and closes the file.
  void close_file();

  // Writes the current state to disk into the opened file.
  void flush();

  // Creates or opens a pre-existing group.
  // Throws if there is a data entry with the same name in the current group.
  void open_group(const std::string& name);
  // Closes the topmost open group.
  // Precondition: the current group is not the root.
  void close_group();

  constexpr static bool is_reader = false;
  constexpr static bool is_writer = true;

  void set_verbose(bool verbose) noexcept {
    verbose_ = verbose;
  }

  // Returns true if a file is open.
  operator bool() const noexcept;

  // Commits to the internal data representation the value of "obj".
  template <class T>
  void execute(const std::string& name, const T& obj) noexcept;
  template <class Scalar, class Domain>
  void execute(const std::string& name, const func::function<Scalar, Domain>& f) noexcept;
  template <class Scalar>
  void execute(const std::string& name, const linalg::Matrix<Scalar, dca::linalg::CPU>& m) noexcept;
  template <class T>
  void execute(const std::string& name, const std::unique_ptr<T>& ptr) noexcept {
    if (ptr)
      execute(name, *ptr);
  }
  template <class T>
  void execute(const T& f) noexcept {
    execute(f.get_name(), f);
  }
  template <class T>
  void execute(const std::unique_ptr<T>& f) noexcept {
    if (f)
      execute(f->get_name(), *f);
  }

private:
  bool verbose_;
  std::stack<details::JSONGroup*> open_groups_;
  details::JSONGroup root_;

  std::ofstream file_;
};

template <class T>
void JSONWriter::execute(const std::string& name, const T& obj) noexcept {
  open_groups_.top()->addEntry(name, obj);
}

template <class Scalar, class Domain>
void JSONWriter::execute(const std::string& name, const func::function<Scalar, Domain>& f) noexcept {
  if (verbose_)
    std::cout << "\t starts writing function : " << f.get_name() << "\n";

  open_group(name);

  execute("name", f.get_name());
  execute("domain-sizes", f.getDomainSizes());
  execute("data", f.getValues());

  close_group();
}

template <class Scalar>
void JSONWriter::execute(const std::string& name,
                         const linalg::Matrix<Scalar, dca::linalg::CPU>& m) noexcept {
  open_group(name);

  std::vector<std::vector<Scalar>> data(m.nrRows());
  for (int i = 0; i < m.nrRows(); ++i) {
    data[i].resize(m.nrCols());
    for (int j = 0; j < m.nrCols(); ++j) {
      data[i][j] = m(i, j);
    }
  }

  execute("name", m.get_name());
  execute("size", m.size());
  execute("data", data);

  close_group();
}

}  // namespace dca::io

#endif  // DCA_IO_JSON_JSON_WRITER_HPP
