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
#ifdef DCA_HAVE_ADIOS2
#include "dca/io/adios2/adios2_writer.hpp"
#endif

#define RETURN_IF_NOT_PARALLEL(x) if(concurrency_.id() != concurrency_.first() && !isADIOS2()) \
                                      return x

namespace dca::io {

template <class Concurrency>
class Writer {
public:
  using CT = Concurrency;
  using DCAWriterVariant = std::variant<io::HDF5Writer, io::JSONWriter
#ifdef DCA_HAVE_ADIOS2
                                        ,
                                        io::ADIOS2Writer<Concurrency>
#endif
                                        >;
  /** Constructor for writer, pretty tortured due to optional Adios2.
   *  \param [in] adios        if build with ADIOS2 support adios2 env reference.
   *  \param [in] concurrency  reference to the concurrency environment
   *  \param [in] format       output format: ADIOS2, HDF5 or JSON.
   *  \param [in] verbose      If true, the writer outputs a short log whenever it is executed.
   *
   *  When using a format other than ADIOS2 the concurrency and adios arguments aren't
   *  used by the writer.
   */
  Writer(
#ifdef DCA_HAVE_ADIOS2
      adios2::ADIOS& adios,
#endif
      Concurrency& concurrency, const std::string& format, bool verbose = true)
      :
#ifdef DCA_HAVE_ADIOS2
        adios_(adios),
#endif
        concurrency_(concurrency) {
    if (format == "HDF5") {
      writer_.template emplace<io::HDF5Writer>(verbose);
    }
    else if (format == "JSON") {
      writer_.template emplace<io::JSONWriter>(verbose);
    }
#ifdef DCA_HAVE_ADIOS2
    else if (format == "ADIOS2") {
      writer_.template emplace<io::ADIOS2Writer<Concurrency>>(adios_, &concurrency_, verbose);
    }
#endif
    else {
      throw(std::logic_error("Invalid output format"));
    }
  }

  constexpr static bool is_reader = false;
  constexpr static bool is_writer = true;

  DCAWriterVariant& getUnderlying() {
    return writer_;
  }

#ifdef DCA_HAVE_ADIOS2
  bool isADIOS2() {
    return std::holds_alternative<io::ADIOS2Writer<Concurrency>>(writer_);
  }
#else
  bool isADIOS2() {
    return false;
  }
#endif

  bool isHDF5() {
    return std::holds_alternative<io::HDF5Writer>(writer_);
  }
  
  bool isOpen() { return is_open_; }
  
  void begin_step() {
    RETURN_IF_NOT_PARALLEL();
    std::visit([&](auto& var) { var.begin_step(); }, writer_);
  }

  void end_step() {
    RETURN_IF_NOT_PARALLEL();
    std::visit([&](auto& var) { var.end_step(); }, writer_);
  }

  void open_file(const std::string& file_name, bool overwrite = true) {
    RETURN_IF_NOT_PARALLEL();
    std::visit([&](auto& var) { var.open_file(file_name, overwrite); }, writer_);
    is_open_ = true;
  }

  void close_file() {
    RETURN_IF_NOT_PARALLEL();
    std::visit([&](auto& var) { var.close_file(); }, writer_);
    is_open_ = false;
  }

  /** For writing open_group is expected to always return true
   */
  bool open_group(const std::string& new_path) {
    RETURN_IF_NOT_PARALLEL(true);
    return std::visit([&](auto& var) -> bool { return var.open_group(new_path); }, writer_);
  }

  void close_group() {
    RETURN_IF_NOT_PARALLEL();
    std::visit([&](auto& var) { var.close_group(); }, writer_);
  }

  void flush() {
    RETURN_IF_NOT_PARALLEL();
    std::visit([&](auto& var) { var.flush(); }, writer_);
  }

  template <class... Args>
  bool execute(const Args&... args) {
    RETURN_IF_NOT_PARALLEL(true);
    return std::visit([&](auto& var) ->bool { return var.execute(args...); }, writer_);
  }

  template <class... Args>
  void executePartial(const Args&... args) {
    RETURN_IF_NOT_PARALLEL();
    std::visit([&](auto& var) { var.executePartial(args...); }, writer_);
  }

  template <class... Args>
  void rewrite(const std::string& name, const Args&... args) {
    if (concurrency_.id() == concurrency_.first()) {
      std::visit(
          [&](auto& var) {
            var.erase(name);
            var.execute(name, args...);
          },
          writer_);
    }
  }

  operator bool() const noexcept {
    return std::visit([&](const auto& var) { return static_cast<bool>(var); }, writer_);
  }

  void lock() {
    RETURN_IF_NOT_PARALLEL();
    mutex_.lock();
  }

  void unlock() {
    RETURN_IF_NOT_PARALLEL();
    mutex_.unlock();
  }

  void set_verbose(bool verbose) {
    std::visit([&](auto& var) { var.set_verbose(verbose); }, writer_);
  }

  Concurrency& get_concurrency() { return concurrency_; }
  
  dca::parallel::thread_traits::mutex_type& get_mutex() { return mutex_; }
private:
  dca::parallel::thread_traits::mutex_type mutex_;
  DCAWriterVariant writer_;
#ifdef DCA_HAVE_ADIOS2
  adios2::ADIOS& adios_;
#endif
  Concurrency& concurrency_;
  bool is_open_;
};

}  // namespace dca::io

#endif  // DCA_IO_WRITER_HPP
