// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Thomas Maier (maierta@ornl.gov)
//
// Template specialization of model parameters for square lattice tow-orbital bilayer Hubbard model for La3Ni2O7.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::twoOrbital<PointGroup>>> {
public:
  ModelParameters()
      : t11_(0.),
        t22_(0.),
        t12_(0.),
        U_(0.),
        V_(0.),
        J_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_t11() const {
    return t11_;
  }
  void set_t11(const double t11) {
    t11_ = t11;
  }

  double get_t22() const {
    return t22_;
  }
  void set_t22(const double t22) {
    t22_ = t22;
  }

  double get_t12() const {
    return t12_;
  }
  void set_t12(const double t12) {
    t12_ = t12;
  }

  double get_U() const {
    return U_;
  }
  void set_U(const double U) {
    U_ = U;
  }

  double get_V() const {
    return V_;
  }
  void set_V(const double V) {
    V_ = V;
  }

  double get_J() const {
    return J_;
  }
  void set_J(const double J) {
    J_ = J;
  }

private:
  double t11_;
  double t22_;
  double t12_;
  double U_;
  double V_;
  double J_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::twoOrbital<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t11_);
  buffer_size += concurrency.get_buffer_size(t22_);
  buffer_size += concurrency.get_buffer_size(t12_);
  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(J_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::twoOrbital<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t11_);
  concurrency.pack(buffer, buffer_size, position, t22_);
  concurrency.pack(buffer, buffer_size, position, t12_);
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, J_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::twoOrbital<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t11_);
  concurrency.unpack(buffer, buffer_size, position, t22_);
  concurrency.unpack(buffer, buffer_size, position, t12_);
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, J_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::twoOrbital<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("two-orbital-model");

    try {
      reader_or_writer.execute("t11", t11_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t22", t22_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t12", t12_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U", U_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("V", V_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("J", J_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}
