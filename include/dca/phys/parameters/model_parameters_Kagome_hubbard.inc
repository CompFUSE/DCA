// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific
// publications.
//
// Author: Peter Doak (doakpw@ornl.govb)
//
// Template specialization of model parameters for Kagome lattice Hubbard model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::KagomeHubbard<PointGroup>>> {
public:
  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_t() const {
    return t_;
  }
  void set_t(double t) {
    t_ = t;
  }

  double get_U() const {
    return U_;
  }
  void set_U(double U) {
    U_ = U;
  }
private:
  double t_ = 0;
  double U_ = 0;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::KagomeHubbard<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_);
  buffer_size += concurrency.get_buffer_size(U_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::KagomeHubbard<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_);
  concurrency.pack(buffer, buffer_size, position, U_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::KagomeHubbard<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_);
  concurrency.unpack(buffer, buffer_size, position, U_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::KagomeHubbard<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("Kagome-Hubbard-model");
  }
  catch (...) {
    throw std::runtime_error("You must specify a Kagome-Hubbard-model seciton in input!");
  }
    try {
      reader_or_writer.execute("t", t_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U", U_);
    }
    catch (const std::exception& r_e) {
    }
    if (std::abs(U_) <= 1e-3)
      throw std::runtime_error("For Kagome model abs(U) must be greater than 1e-3");

    reader_or_writer.close_group();
}

