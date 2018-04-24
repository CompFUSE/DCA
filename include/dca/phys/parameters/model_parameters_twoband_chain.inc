// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization of model parameters for square lattice bilayer Hubbard model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::twoband_chain<PointGroup>>> {
public:
  ModelParameters() : tx_(0.), ty_(0.), t_perp_(1.), U_(0.), V_(0.), V_prime_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  auto get_t() const {
    return std::array<double, 2>{tx_, ty_};
  }
  void set_t(const double t) {
    tx_ = t;
  }

  double get_t_prime() const {
    return ty_;
  }
  void set_t_prime(const double t_prime) {
    ty_ = t_prime;
  }

  double get_t_perp() const {
    return t_perp_;
  }
  void set_t_perp(const double t_perp) {
    t_perp_ = t_perp;
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

  double get_V_prime() const {
    return V_prime_;
  }
  void set_V_prime(const double V_prime) {
    V_prime_ = V_prime;
  }

private:
  double tx_;
  double ty_;
  double t_perp_;
  double U_;
  double V_;
  double V_prime_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::twoband_chain<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(tx_);
  buffer_size += concurrency.get_buffer_size(ty_);
  buffer_size += concurrency.get_buffer_size(t_perp_);
  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::twoband_chain<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, tx_);
  concurrency.pack(buffer, buffer_size, position, ty_);
  concurrency.pack(buffer, buffer_size, position, t_perp_);
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::twoband_chain<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, tx_);
  concurrency.unpack(buffer, buffer_size, position, ty_);
  concurrency.unpack(buffer, buffer_size, position, t_perp_);
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::twoband_chain<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("twoband-chain-model");
    try {
      reader_or_writer.execute("t", tx_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("tx", tx_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ty", ty_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t-perp", t_perp_);
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
      reader_or_writer.execute("V-prime", V_prime_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}
