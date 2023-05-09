// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// Wrapper to an instance of ADIOS2Reader, HDF5Reader or JSONReader.

#ifndef DCA_IO_READER_HPP
#define DCA_IO_READER_HPP

#include <mutex>
#include <string>
#include <variant>

#include "dca/config/haves_defines.hpp"
#include "dca/platform/dca_gpu.h"

#include "dca/io/io_types.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/json/json_reader.hpp"

#ifdef DCA_HAVE_ADIOS2
#include "dca/io/adios2/adios2_reader.hpp"
#endif

namespace dca::io {

template <class Concurrency>
class Reader {
public:
  using DCAReaderVariant = std::variant<io::HDF5Reader, io::JSONReader
#ifdef DCA_HAVE_ADIOS2
                                        ,
                                        io::ADIOS2Reader<Concurrency>
#endif
                                        >;
  /**
   * \param[in] concureency   reference to current concurrency env
   * \param[in] format        format as IOType
   * \param[in] verbose       If true, reader does some logging
   */
  Reader(const Concurrency& concurrency, const IOType format, bool verbose = true)
      : concurrency_(concurrency) {
    switch (format) {
      case IOType::HDF5:
        reader_.template emplace<io::HDF5Reader>(verbose);
        break;
      case IOType::JSON:
        reader_.template emplace<io::JSONReader>(verbose);
        break;
#ifdef DCA_HAVE_ADIOS2
      case IOType::ADIOS2:
        reader_.template emplace<io::ADIOS2Reader<Concurrency>>(&concurrency, verbose);
        break;
#endif
    }
  }

  /** DEPRECATED -- Support for format type as string 
   * \param[in] format        string representation of format, since parameters still use string
   *
   * \todo remove need for this constructor store IO format as enum class type.
   */
  Reader(const Concurrency& concurrency, const std::string& format, bool verbose = true)
      : Reader(concurrency, stringToIOType(format), verbose) {}

#ifdef DCA_HAVE_ADIOS2
  Reader(adios2::ADIOS& adios, const Concurrency& concurrency, const std::string& format,
         bool verbose = true)
      : concurrency_(concurrency) {
    if (format == "HDF5") {
      throw(std::logic_error("ADIOS2 reference not an argument for hdf5 reader"));
    }
    else if (format == "JSON") {
      throw(std::logic_error("ADIOS2 reference not an argument for json reader"));
    }
    else if (format == "ADIOS2") {
      reader_.template emplace<io::ADIOS2Reader<Concurrency>>(adios, &concurrency, verbose);
    }
    else {
      throw(std::logic_error("Invalid input format"));
    }
  }
#endif

  constexpr static bool is_reader = true;
  constexpr static bool is_writer = false;

  void open_file(const std::string& file_name) {
    std::visit([&](auto& var) { var.open_file(file_name); }, reader_);
  }
  void close_file() {
    std::visit([&](auto& var) { var.close_file(); }, reader_);
  }

  /** For reading input there is great utility in knowing if a group is present.
   *  It isn't an exceptional circumstance if a group is not present.
   */
  bool open_group(const std::string& new_path) {
    return std::visit([&](auto& var) -> bool { return var.open_group(new_path); }, reader_);
  }

  void close_group() {
    std::visit([&](auto& var) { var.close_group(); }, reader_);
  }

  std::size_t getStepCount() {
    return std::visit([&](auto& var) ->std::size_t { return var.getStepCount(); }, reader_);
  }
  
  void begin_step() {
    std::visit([&](auto& var) { var.begin_step(); }, reader_);
  }

  void end_step() {
    std::visit([&](auto& var) { var.end_step(); }, reader_);
  }

  std::string get_path() {
    return std::visit([&](auto& var) -> std::string { return var.get_path(); }, reader_);
  }
  
  template <class... Args>
  bool execute(Args&&... args) noexcept {
    return std::visit([&](auto& var) -> bool { return var.execute(std::forward<Args>(args)...); },
                      reader_);
  }

  DCAReaderVariant& getUnderlying() {
    return reader_;
  }

private:
  DCAReaderVariant reader_;
  const Concurrency& concurrency_;
};

extern template class Reader<dca::parallel::NoConcurrency>;
#ifdef DCA_HAVE_MPI
extern template class Reader<dca::parallel::MPIConcurrency>;
#endif

}  // namespace dca::io

#endif  // DCA_IO_READER_HPP
