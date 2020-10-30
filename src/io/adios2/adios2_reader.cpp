// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Norbert podhorszki (pnorbert@ornl.gov)
//
// This file implements adios2_reader.hpp.

#include "dca/io/adios2/adios2_reader.hpp"

#include <fstream>
#include <stdexcept>

namespace dca {
namespace io {
// dca::io::

template <class CT>
ADIOS2Reader<CT>::ADIOS2Reader(adios2::ADIOS& adios, const CT* concurrency,
                           bool verbose)
    : adios_(adios),
      verbose_(verbose),
      concurrency_(concurrency) {}

template <class CT>
ADIOS2Reader<CT>::~ADIOS2Reader() {
  if (file_)
    close_file();
}


template <class CT>
void ADIOS2Reader<CT>::open_file(const std::string& file_name) {
  if (verbose_) {
    std::cout << "\t ADIOS2Reader: Open for Read file : " << file_name << "\n";
  }

  io_name_ = file_name;
  file_name_ = file_name;
  io_ = adios_.DeclareIO(io_name_);
  file_ = io_.Open(file_name_, adios2::Mode::Read);
}

template <class CT>
void ADIOS2Reader<CT>::close_file() {
  if (file_) {
    file_.Close();
    adios_.RemoveIO(io_name_);
  }
}

template <class CT>
void ADIOS2Reader<CT>::open_group(const std::string& name) {
  size_t len = name.size();
  // remove trailing / from name
  for (; name[len - 1] == '/'; --len)
    ;

  my_paths_.push_back(std::string(name, 0, len));
}

template <class CT>
void ADIOS2Reader<CT>::close_group() {
  my_paths_.pop_back();
}

template <class CT>
std::string ADIOS2Reader<CT>::get_path(const std::string& name) {
  std::string path = "/";

  for (size_t i = 0; i < my_paths_.size(); i++) {
    path += my_paths_[i] + "/";
  }

  if (!name.empty()) {
    path += name;
  }

  return path;
}

template <class CT>
bool ADIOS2Reader<CT>::execute(const std::string& name, std::string& value) {
  std::string full_name = get_path(name);
  if (!exists(full_name)) {
    return false;
  }

  /* ADIOS will resize string to match size of incoming data */
  file_.Get<std::string>(full_name, value, adios2::Mode::Sync);

  return true;
}

template <class CT>
bool ADIOS2Reader<CT>::execute(const std::string& name, std::vector<std::string>& value) {
  std::string full_name = get_path(name);
  bool retval = true;
  if (exists(full_name)) {
    value.resize(1);
    /* ADIOS will resize string to match size of incoming data */
    file_.Get<std::string>(full_name, value[0], adios2::Mode::Sync);
  }
  else {
    auto strvecAttr = io_.InquireAttribute<std::string>(full_name);
    if (strvecAttr) {
      const std::vector<std::string>& strings = strvecAttr.Data();
      value.resize(strings.size());
      for (int i = 0; i < strings.size(); ++i) {
        value[i] = strings[i];
      }
    }
    else {
      retval = false;
    }
  }

  return retval;
}

template <class CT>
bool ADIOS2Reader<CT>::exists(const std::string& name) const {
  std::string varType = io_.VariableType(name);
  return !varType.empty();
}

template class ADIOS2Reader<dca::parallel::NoConcurrency>;
#ifdef DCA_HAVE_MPI
template class ADIOS2Reader<dca::parallel::MPIConcurrency>;
#endif

}  // namespace io
}  // namespace dca
