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

#include "dca/util/type_help.hpp"
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
  using DCAReaderVariant = std::variant<std::monostate,
#ifdef DCA_HAVE_ADIOS2
                                        io::ADIOS2Reader<Concurrency>,
#endif
                                        io::HDF5Reader, io::JSONReader>;
  /**
   * \param[in] concureency   reference to current concurrency env
   * \param[in] format        format as IOType
   * \param[in] verbose       If true, reader does some logging
   */
  Reader(Concurrency& concurrency, const std::string& format_tag, bool verbose = true)
    : Reader(concurrency, stringToIOType(format_tag), verbose) {}
  
  Reader(Concurrency& concurrency, IOType format, bool verbose = true)
      : concurrency_(concurrency) {
    switch (format) {
      case IOType::HDF5:
        reader_.template emplace<io::HDF5Reader>(verbose);
        break;
      case IOType::JSON:
        reader_.template emplace<io::JSONReader>(verbose);
        break;
      case IOType::ADIOS2:
#ifdef DCA_HAVE_ADIOS2
        reader_.template emplace<io::ADIOS2Reader<Concurrency>>(concurrency, verbose);

#endif
        break;
    }
  }

  
  constexpr static bool is_reader = true;
  constexpr static bool is_writer = false;

  void open_file(const std::string& file_name) {
    std::visit(
        [&](auto& var) {
          if constexpr (std::is_same_v<std::monostate&, decltype(var)>)
            throw std::runtime_error(
                "No operations should ever occur on monostate in reader variant");
          else
            var.open_file(file_name);
        },
        reader_);
  }
  void close_file() {
    std::visit(
        [&](auto& var) {
          if constexpr (std::is_same_v<std::monostate&, decltype(var)>)
            throw std::runtime_error(
                "No operations should ever occur on monostate in reader variant");
          else
            var.close_file();
        },
        reader_);
  }

  /** For reading input there is great utility in knowing if a group is present.
   *  It isn't an exceptional circumstance if a group is not present.
   */
  bool open_group(const std::string& new_path) {
    return std::visit(
        [&](auto& var) -> bool {
          if constexpr (std::is_same_v<std::monostate&, decltype(var)>)
            throw std::runtime_error(
                "No operations should ever occur on monostate in reader variant");
          else
            return var.open_group(new_path);
        },
        reader_);
  }

  void close_group() {
    std::visit(
        [&](auto& var) {
          if constexpr (std::is_same_v<std::monostate&, decltype(var)>)
            throw std::runtime_error(
                "No operations should ever occur on monostate in reader variant");
          else
            var.close_group();
        },
        reader_);
  }

  long getStepCount() {
    return std::visit(
        [&](auto& var) -> std::size_t {
          if constexpr (std::is_same_v<std::monostate&, decltype(var)>)
            throw std::runtime_error(
                "No operations should ever occur on monostate in reader variant");
          else
            return var.getStepCount();
        },
        reader_);
  }

  void begin_step() {
    std::visit(
        [&](auto& var) {
          if constexpr (std::is_same_v<std::monostate&, decltype(var)>)
            throw std::runtime_error(
                "No operations should ever occur on monostate in reader variant");
          else
            var.begin_step();
        },
        reader_);
  }

  void end_step() {
    std::visit(
        [&](auto& var) {
          if constexpr (std::is_same_v<std::monostate&, decltype(var)>)
            throw std::runtime_error(
                "No operations should ever occur on monostate in reader variant");
          else
            var.end_step();
        },
        reader_);
  }

  std::string get_path() {
    return std::visit(
        [&](auto& var) -> std::string {
          if constexpr (std::is_same_v<std::monostate&, decltype(var)>)
            throw std::runtime_error(
                "No operations should ever occur on monostate in reader variant");
          else
            return var.get_path();
        },
        reader_);
  }

  template <class... Args>
  bool execute(Args&&... args) noexcept {
    return std::visit(
        [&](auto& var) -> bool {
          if constexpr (std::is_same_v<std::monostate&, decltype(var)>)
            throw std::runtime_error(
                "No operations should ever occur on monostate in reader variant");
          else
            return var.execute(std::forward<Args>(args)...);
        },
        reader_);
  }

  DCAReaderVariant& getUnderlying() {
    return reader_;
  }

private:
  DCAReaderVariant reader_;
  Concurrency& concurrency_;
};

extern template class Reader<dca::parallel::NoConcurrency>;
#ifdef DCA_HAVE_MPI
extern template class Reader<dca::parallel::MPIConcurrency>;
#endif

}  // namespace dca::io

#endif  // DCA_IO_READER_HPP
