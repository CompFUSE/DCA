// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Andrei Plamada (plamada@itp.phys.ethz.ch)
//
// Template specialization of model parameters for 2D 2-band model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::twoband_lattice<PointGroup>>> {
public:
  ModelParameters() : ei0_(0.), eb0_(0.), t0_(0.), U0_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_ei0() const {
    return ei0_;
  }
  double get_eb0() const {
    return eb0_;
  }
  double get_t0() const {
    return t0_;
  }
  double get_U0() const {
    return U0_;
  }

private:
  double ei0_;
  double eb0_;
  double t0_;
  double U0_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::twoband_lattice<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(ei0_);
  buffer_size += concurrency.get_buffer_size(eb0_);
  buffer_size += concurrency.get_buffer_size(t0_);
  buffer_size += concurrency.get_buffer_size(U0_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::twoband_lattice<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, ei0_);
  concurrency.pack(buffer, buffer_size, position, eb0_);
  concurrency.pack(buffer, buffer_size, position, t0_);
  concurrency.pack(buffer, buffer_size, position, U0_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::twoband_lattice<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, ei0_);
  concurrency.unpack(buffer, buffer_size, position, eb0_);
  concurrency.unpack(buffer, buffer_size, position, t0_);
  concurrency.unpack(buffer, buffer_size, position, U0_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::twoband_lattice<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("twoband-model");

    try {
      reader_or_writer.execute("ei0", ei0_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("eb0", eb0_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t0", t0_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U0", U0_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\nNo 2-band model parameters defined!" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  std::cout << "\n\n"
            << "\t twoband-model : \n"
            << "\t--------------- \n\n"
            << "\t\t ei0  : " << ei0_ << "\n"
            << "\t\t eb0  : " << eb0_ << "\n"
            << "\t\t  t0  : " << t0_ << "\n"
            << "\t\t  U0  : " << U0_ << "\n"
            << "\n"
            << std::endl;
}
