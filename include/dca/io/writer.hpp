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

#include <string>
#include <variant>

#include "dca/config/threading.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_writer.hpp"
#ifdef DCA_WITH_ADIOS2
#include "dca/io/adios2/adios2_writer.hpp"
#endif

namespace dca::io {

template <class Concurrency>
class Writer {
public:
  // In: format. output format, HDF5 or JSON.
  // In: verbose. If true, the writer outputs a short log whenever it is executed.
  Writer(Concurrency& concurrency, const std::string& format, bool verbose = true)
      : concurrency_(concurrency) {
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
    // currently only the ADIOS2Writer supports parallel writes
#ifdef DCA_HAVE_ADIOS2
    if constexpr (std::is_same<decltype(writer_), ADIOS2Writer<Concurrency>>::value) {
      std::visit([&](auto& var) { var.execute(args...); }, writer_);
    }
    else {
#endif
      if (concurrency_.id() == concurrency_.first()) {
        std::visit([&](auto& var) { var.execute(args...); }, writer_);
      }
#ifdef DCA_HAVE_ADIOS2
    }
#endif
  }

  template <class... Args>
  void rewrite(const std::string& name, const Args&... args) {
#ifdef DCA_HAVE_ADIOS2
    if constexpr (std::is_same<decltype(writer_), ADIOS2Writer<Concurrency>>::value) {
      std::visit(
          [&](auto& var) {
            var.erase(name);
            var.execute(name, args...);
          },
          writer_);
    }
    else {
#endif
      if (concurrency_.id() == concurrency_.first()) {
        std::visit(
            [&](auto& var) {
              var.erase(name);
              var.execute(name, args...);
            },
            writer_);
      }
#ifdef DCA_HAVE_ADIOS2
    }
#endif
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
  dca::parallel::thread_traits::mutex_type mutex_;
  std::variant<io::HDF5Writer, io::JSONWriter
#ifdef DCA_HAVE_ADIOS2
               ,io::ADIOS2Writer
#endif
               > writer_;
  Concurrency& concurrency_;
};

}  // namespace dca::io

#endif  // DCA_IO_WRITER_HPP
