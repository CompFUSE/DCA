// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific
// publications.
//
// Author: Thomas Maier (maierta@ornl.govb)
//
// Template specialization of model parameters for two-orbital altermagnet Hubbard model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::Altermagnet<PointGroup>>> {
public:
  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_t0() const {
    return t0_;
  }
  void set_t0(double t0) {
    t0_ = t0;
  }

  double get_t1() const {
    return t1_;
  }
  void set_t1(double t1) {
    t0_ = t1;
  }

  double get_t2() const {
    return t2_;
  }
  void set_t2(double t2) {
    t0_ = t2;
  }


  double get_U() const {
    return U_;
  }
  void set_U(double U) {
    U_ = U;
  }
private:
  double t0_ = 0;
  double t1_ = 0;
  double t2_ = 0;
  double U_ = 0;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::Altermagnet<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t0_);
  buffer_size += concurrency.get_buffer_size(t1_);
  buffer_size += concurrency.get_buffer_size(t2_);
  buffer_size += concurrency.get_buffer_size(U_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::Altermagnet<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t0_);
  concurrency.pack(buffer, buffer_size, position, t1_);
  concurrency.pack(buffer, buffer_size, position, t2_);
  concurrency.pack(buffer, buffer_size, position, U_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::Altermagnet<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t0_);
  concurrency.unpack(buffer, buffer_size, position, t1_);
  concurrency.unpack(buffer, buffer_size, position, t2_);
  concurrency.unpack(buffer, buffer_size, position, U_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::Altermagnet<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("Altermagnet-Hubbard-model");
  }
  catch (...) {
    throw std::runtime_error("You must specify a Altermagnet-Hubbard-model section in input!");
  }
    try {
      reader_or_writer.execute("t0", t0_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t1", t1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t2", t2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U", U_);
    }
    catch (const std::exception& r_e) {
    }
    if (std::abs(U_) <= 1e-3)
      throw std::runtime_error("For Altermagnet model abs(U) must be greater than 1e-3");

    reader_or_writer.close_group();
}

