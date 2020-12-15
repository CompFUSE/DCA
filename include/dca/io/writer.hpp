// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Wrapper to an instance of HDF5Writer of JSONWriter.

#ifndef DCA_IO_WRITER_HPP
#define DCA_IO_WRITER_HPP

#include <mutex>
#include <string>
#include <variant>

#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_writer.hpp"

namespace dca::io {

class Writer {
public:
  // In: format. output format, HDF5 or JSON.
  // In: verbose. If true, the writer outputs a short log whenever it is executed.
  Writer(const std::string& format, bool verbose = true) {
    if (format == "HDF5") {
      writer_.emplace<io::HDF5Writer>(verbose);
    }
    else if (format == "JSON") {
      writer_.emplace<io::JSONWriter>(verbose);
    }
    else {
      throw(std::logic_error("Invalid output format"));
    }
  }

  constexpr static bool is_reader = false;
  constexpr static bool is_writer = true;

  void open_file(const std::string& file_name, bool overwrite = true) {
    std::visit([&](auto& var) { var.open_file(file_name, overwrite); }, writer_);
  }

  void close_file() {
    std::visit([&](auto& var) { var.close_file(); }, writer_);
  }

  void open_group(const std::string& new_path) {
    std::visit([&](auto& var) { var.open_group(new_path); }, writer_);
  }
  void close_group() {
    std::visit([&](auto& var) { var.close_group(); }, writer_);
  }

  template <class... Args>
  void execute(const Args&... args) {
    std::visit([&](auto& var) { var.execute(args...); }, writer_);
  }

  operator bool() const noexcept {
    return std::visit([&](const auto& var) { return static_cast<bool>(var); }, writer_);
  }

  void lock() {
    mutex_.lock();
  }

  void unlock() {
    mutex_.unlock();
  }

  void set_verbose(bool verbose) {
    std::visit([&](auto& var) { var.set_verbose(verbose); }, writer_);
  }

private:
  std::mutex mutex_;
  std::variant<io::HDF5Writer, io::JSONWriter> writer_;
};

}  // namespace dca::io

#endif  // DCA_IO_WRITER_HPP
