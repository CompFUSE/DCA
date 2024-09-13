// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization of model parameters for square lattice bilayer Hubbard model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>> {
public:
  ModelParameters() : t_(0.), t_prime_(0.), t_perp_(0.), U_(0.), V_(0.), V_prime_(0.) {}

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
  void set_t(const double t) {
    t_ = t;
  }

  double get_t_prime() const {
    return t_prime_;
  }
  void set_t_prime(const double t_prime) {
    t_prime_ = t_prime;
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
  double t_;
  double t_prime_;
  double t_perp_;
  double U_;
  double V_;
  double V_prime_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_);
  buffer_size += concurrency.get_buffer_size(t_prime_);
  buffer_size += concurrency.get_buffer_size(t_perp_);
  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_);
  concurrency.pack(buffer, buffer_size, position, t_prime_);
  concurrency.pack(buffer, buffer_size, position, t_perp_);
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_);
  concurrency.unpack(buffer, buffer_size, position, t_prime_);
  concurrency.unpack(buffer, buffer_size, position, t_perp_);
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("bilayer-Hubbard-model");

    try {
      reader_or_writer.execute("t", t_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t-prime", t_prime_);
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
    std::cout << "\n  bilayer-Hubbard-model parameter group is not defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  if (ReaderOrWriter::is_reader) {
    std::cout << "\n\n"
              << "\tbilayer-Hubbard-model : \n"
              << "\t--------------- \n\n"
              << "\t\t t  : " << t_ << "\n"
              << "\t\t t_prime  : " << t_prime_ << "\n"
              << "\t\t t_perp  : " << t_perp_ << "\n"
              << "\t\t U : " << U_ << "\n"
              << "\t\t V : " << V_ << "\n"
              << "\t\t V-prime : " << V_prime_ << "\n"
              << "\n"
              << std::endl;
  }

}
