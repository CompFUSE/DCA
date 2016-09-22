// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class reads, stores, and writes the exact diagonalization solver parameters.

#ifndef DCA_PHYS_PARAMETERS_ED_SOLVER_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_ED_SOLVER_PARAMETERS_HPP

#include <stdexcept>
#include <string>

namespace dca {
namespace phys {
namespace params {
// dca::phys::params::

class EdSolverParameters {
public:
  EdSolverParameters()
      : eigenvalue_cut_off_(1.e-6),
        ed_method_("default"),
        occupation_(0),
        magnetization_(0),
        check_orthogonality_of_states_("false") {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_eigenvalue_cut_off() const {
    return eigenvalue_cut_off_;
  }
  const std::string& get_ed_method() const {
    return ed_method_;
  }
  int get_occupation() const {
    return occupation_;
  }
  int get_magnetization() const {
    return magnetization_;
  }
  bool check_orthogonality_of_states() const {
    return (check_orthogonality_of_states_ == "true");
  }

private:
  double eigenvalue_cut_off_;
  std::string ed_method_;
  int occupation_;
  int magnetization_;
  std::string check_orthogonality_of_states_;
};

template <typename Concurrency>
int EdSolverParameters::getBufferSize(const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(eigenvalue_cut_off_);
  buffer_size += concurrency.get_buffer_size(ed_method_);
  buffer_size += concurrency.get_buffer_size(occupation_);
  buffer_size += concurrency.get_buffer_size(magnetization_);
  buffer_size += concurrency.get_buffer_size(check_orthogonality_of_states_);

  return buffer_size;
}

template <typename Concurrency>
void EdSolverParameters::pack(const Concurrency& concurrency, int* buffer, int buffer_size,
                              int& position) const {
  concurrency.pack(buffer, buffer_size, position, eigenvalue_cut_off_);
  concurrency.pack(buffer, buffer_size, position, ed_method_);
  concurrency.pack(buffer, buffer_size, position, occupation_);
  concurrency.pack(buffer, buffer_size, position, magnetization_);
  concurrency.pack(buffer, buffer_size, position, check_orthogonality_of_states_);
}

template <typename Concurrency>
void EdSolverParameters::unpack(const Concurrency& concurrency, int* buffer, int buffer_size,
                                int& position) {
  concurrency.unpack(buffer, buffer_size, position, eigenvalue_cut_off_);
  concurrency.unpack(buffer, buffer_size, position, ed_method_);
  concurrency.unpack(buffer, buffer_size, position, occupation_);
  concurrency.unpack(buffer, buffer_size, position, magnetization_);
  concurrency.unpack(buffer, buffer_size, position, check_orthogonality_of_states_);
}

template <typename ReaderOrWriter>
void EdSolverParameters::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("ED-solver-parameters");

    try {
      reader_or_writer.execute("eigenvalue-cut-off", eigenvalue_cut_off_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ED-method", ed_method_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("occupation", occupation_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("magnetization", magnetization_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("check_orthogonality_of_states", check_orthogonality_of_states_);
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
