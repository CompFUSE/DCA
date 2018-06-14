// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the exact diagonalization solver parameters.

#ifndef DCA_PHYS_PARAMETERS_ED_SOLVER_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_ED_SOLVER_PARAMETERS_HPP

#include <stdexcept>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class EdSolverParameters {
public:
  EdSolverParameters() : eigenvalue_cut_off_(1.e-6) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_eigenvalue_cut_off() const {
    return eigenvalue_cut_off_;
  }

private:
  double eigenvalue_cut_off_;
};

template <typename Concurrency>
int EdSolverParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;
  buffer_size += concurrency.get_buffer_size(eigenvalue_cut_off_);
  return buffer_size;
}

template <typename Concurrency>
void EdSolverParameters::pack(const Concurrency& concurrency, char* buffer, int buffer_size,
                              int& position) const {
  concurrency.pack(buffer, buffer_size, position, eigenvalue_cut_off_);
}

template <typename Concurrency>
void EdSolverParameters::unpack(const Concurrency& concurrency, char* buffer, int buffer_size,
                                int& position) {
  concurrency.unpack(buffer, buffer_size, position, eigenvalue_cut_off_);
}

template <typename ReaderOrWriter>
void EdSolverParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("ED");

    try {
      reader_or_writer.execute("eigenvalue-cut-off", eigenvalue_cut_off_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}

}  // params
}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_ED_SOLVER_PARAMETERS_HPP
