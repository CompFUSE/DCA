
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Norbert Podhorszki (pnorbert@ornl.gov)
//
// This file implements adios2_writer.hpp.

#include "dca/io/adios2/adios2_writer.hpp"

#include <iostream>
#include <fstream>
#include <stdexcept>

namespace dca {
namespace io {
// dca::io::


template <class Concurrency>
ADIOS2Writer<Concurrency>::ADIOS2Writer(adios2::ADIOS& adios, const Concurrency* concurrency,
                                        bool verbose)
    : adios_(adios), verbose_(verbose), concurrency_(concurrency) {
  if constexpr (std::is_same<decltype(concurrency_), dca::parallel::MPIConcurrency>::value) {
    std::cout << "AdiosWriter on MPI world:" << concurrency_->get_world_id()
       << " of size: " << concurrency_->get_world_size() << '\n';
    std::cout << "comm size: " << concurrency_->size() << " id: " << concurrency_->id() << '\n';
  }
  
}

template <class Concurrency>
ADIOS2Writer<Concurrency>::~ADIOS2Writer() {
  if (file_)
    close_file();
}

template <class Concurrency>
void ADIOS2Writer<Concurrency>::open_file(const std::string& file_name_ref, bool overwrite) {
  adios2::Mode mode = (overwrite ? adios2::Mode::Write : adios2::Mode::Append);
  if (verbose_) {
    std::cout << "\t ADIOS2Writer: Open for " << (overwrite ? "Write" : "Append")
              << " file : " << file_name_ref << "\n";
    std::cout << "On world";
  }

  io_name_ = file_name_ref;
  file_name_ = file_name_ref;
  io_ = adios_.DeclareIO(io_name_);
  file_ = io_.Open(file_name_, mode);
  // This is true if m_isClosed is false, that doesn't mean the "file" is open.
  if (!file_) {
    std::ostringstream error_message;
    error_message << "ADIOS2Writer<Concurrency>::open_file failed to open " << file_name_ref;
    throw std::ios_base::failure(error_message.str());
  }
}

template <class Concurrency>
void ADIOS2Writer<Concurrency>::close_file() {
  if (file_) {
    file_.Close();
    //This combined with overwrite seems to create issues.
    adios_.RemoveIO(io_name_);
  }
}

template <class Concurrency>
void ADIOS2Writer<Concurrency>::open_group(const std::string& name) {
  size_t len = name.size();
  // remove trailing / from name
  for (; name[len - 1] == '/'; --len)
    ;
  my_paths_.push_back(std::string(name, 0, len));
}

template <class Concurrency>
void ADIOS2Writer<Concurrency>::close_group() {
  my_paths_.pop_back();
}

template <class Concurrency>
std::string ADIOS2Writer<Concurrency>::get_path(const std::string& name) {
  std::string path = "/";

  for (size_t i = 0; i < my_paths_.size(); i++) {
    path += my_paths_[i] + "/";
  }

  if (!name.empty()) {
    path += name;
  }

  return path;
}

template <class Concurrency>
void ADIOS2Writer<Concurrency>::execute(const std::string& name, const std::string& value) {
  std::string full_name = get_path(name);
  if (value.size() == 0) {
    write(full_name, "");
  }
  else {
    write(full_name, value);
  }
}

template <class Concurrency>
void ADIOS2Writer<Concurrency>::execute(const std::string& name,
                                        const std::vector<std::string>& value) {
  std::string full_name = get_path(name);
  if (value.size() == 0) {
    write(full_name, "");
  }
  else if (value.size() == 1) {
    write(full_name, value[0]);
  }
  else {
    // Store vector of string as attribute since it is not supported as variable
    io_.DefineAttribute(full_name, value.data(), value.size());
  }
}

template <class Concurrency>
void ADIOS2Writer<Concurrency>::write(const std::string& name, const std::string& data) {
  adios2::Variable<std::string> v;
  v = io_.DefineVariable<std::string>(name);
  file_.Put(name, data);
}

template <class Concurrency>
void ADIOS2Writer<Concurrency>::addAttribute(const std::string& set, const std::string& name,
                                             const std::string& value) {
  io_.DefineAttribute(name, value, set);
}

template class ADIOS2Writer<dca::parallel::NoConcurrency>;
template class ADIOS2Writer<dca::parallel::MPIConcurrency>;

}  // namespace io
}  // namespace dca
