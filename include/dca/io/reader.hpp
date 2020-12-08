// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Wrapper to an instance of HDF5Reader of JSONReader.

#ifndef DCA_IO_READER_HPP
#define DCA_IO_READER_HPP

#include <mutex>
#include <string>
#include <variant>

#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/json/json_reader.hpp"

namespace dca::io {

class Reader {
public:
  // In: format. output format, HDF5 or JSON.
  // In: verbose. If true, the reader outputs a short log whenever it is executed.
  Reader(const std::string& format, bool verbose = true) {
    if (format == "HDF5") {
      reader_.emplace<io::HDF5Reader>(verbose);
    }
    else if (format == "JSON") {
      reader_.emplace<io::JSONReader>(verbose);
    }
    else {
      throw(std::logic_error("Invalid input format"));
    }
  }

  constexpr static bool is_reader = true;
  constexpr static bool is_writer = false;

  void open_file(const std::string& file_name) {
    std::visit([&](auto& var) { var.open_file(file_name); }, reader_);
  }

  void close_file() {
    std::visit([&](auto& var) { var.close_file(); }, reader_);
  }

  void open_group(const std::string& new_path) {
    std::visit([&](auto& var) { var.open_group(new_path); }, reader_);
  }
  void close_group() {
    std::visit([&](auto& var) { var.close_group(); }, reader_);
  }

  template <class... Args>
  bool execute(Args&&... args) noexcept {
    return std::visit([&](auto& var) -> bool { return var.execute(std::forward<Args>(args)...); },
                      reader_);
  }

private:
  std::variant<io::HDF5Reader, io::JSONReader> reader_;
};

}  // namespace dca::io

#endif  // DCA_IO_READER_HPP
