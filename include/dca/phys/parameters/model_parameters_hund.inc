// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Template specialization of model parameters for a 2D bilayer model with Hund's interaction.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::HundLattice<PointGroup>>> {
public:
  ModelParameters() : t_(1.), t_prime_(0.), t_perp_(0.), U_(4.), V_(0.), V_prime_(0.) {}

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

  double get_Jh() const {
    return Jh_;
  }
  void set_Jh(const double Jh) {
    Jh_ = Jh;
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
  double t_ = 1.;
  double t_prime_ = 0.;
  double t_perp_ = 0.;
  double U_ = 0.;
  double V_ = 0.;
  double V_prime_ = 0.;
  double Jh_ = 0.;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::HundLattice<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_);
  buffer_size += concurrency.get_buffer_size(t_prime_);
  buffer_size += concurrency.get_buffer_size(t_perp_);
  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);
  buffer_size += concurrency.get_buffer_size(Jh_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::HundLattice<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_);
  concurrency.pack(buffer, buffer_size, position, t_prime_);
  concurrency.pack(buffer, buffer_size, position, t_perp_);
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
  concurrency.pack(buffer, buffer_size, position, Jh_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::HundLattice<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_);
  concurrency.unpack(buffer, buffer_size, position, t_prime_);
  concurrency.unpack(buffer, buffer_size, position, t_perp_);
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
  concurrency.unpack(buffer, buffer_size, position, Jh_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::HundLattice<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  auto tryToRead = [&](const std::string& name, auto& obj) {
    try {
      reader_or_writer.execute(name, obj);
    }
    catch (const std::exception& /*err*/) {
    }
  };

  try {
    reader_or_writer.open_group("Hund-model");
    tryToRead("t", t_);
    tryToRead("t-prime", t_prime_);
    tryToRead("t-perp", t_perp_);
    tryToRead("U", U_);
    tryToRead("V", V_);
    tryToRead("V-prime", V_prime_);
    tryToRead("Jh", Jh_);
    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\n  Hund-model parameter group is not defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  if constexpr (ReaderOrWriter::is_reader)
    std::cout << "\n\n"
              << "\tHund-model : \n"
              << "\t--------------- \n\n"
              << "\t\t t  : " << t_ << "\n"
              << "\t\t t' : " << t_prime_ << "\n"
              << "\t\t tz : " << t_perp_ << "\n"
              << "\t\t U  : " << U_ << "\n"
              << "\t\t V  : " << V_ << "\n"
              << "\t\t V' : " << V_prime_ << "\n"
              << "\t\t Jh' : " << Jh_ << "\n"
              << "\n"
              << std::endl;
}
