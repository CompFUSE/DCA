// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization of model parameters for 2D bilayer model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>> {
public:
  ModelParameters()
      : t_(1.), t_prime_(0.), t_perp_(0.), U_(4.), U_prime_(0.), V_(0.), V_prime_(0.) {}

  template <typename Concurrency>
  int getBufferSize(Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_t() const {
    return t_;
  }
  double get_t_prime() const {
    return t_prime_;
  }
  double get_t_perp() const {
    return t_perp_;
  }
  double get_U() const {
    return U_;
  }
  double get_U_prime() const {
    return U_prime_;
  }
  double get_V() const {
    return V_;
  }
  double get_V_prime() const {
    return V_prime_;
  }

private:
  double t_;
  double t_prime_;
  double t_perp_;
  double U_;
  double U_prime_;
  double V_;
  double V_prime_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::getBufferSize(
    Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_);
  buffer_size += concurrency.get_buffer_size(t_prime_);
  buffer_size += concurrency.get_buffer_size(t_perp_);
  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(U_prime_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::pack(
    Concurrency& concurrency, int* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_);
  concurrency.pack(buffer, buffer_size, position, t_prime_);
  concurrency.pack(buffer, buffer_size, position, t_perp_);
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, U_prime_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::unpack(
    Concurrency& concurrency, int* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_);
  concurrency.unpack(buffer, buffer_size, position, t_prime_);
  concurrency.unpack(buffer, buffer_size, position, t_perp_);
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, U_prime_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("bilayer-model");

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
      reader_or_writer.execute("tz", t_perp_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U", U_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U-prime", U_prime_);
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
    std::cout << "\n No bilayer-model parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  std::cout << "\n\n"
            << "\tbilayer-model : \n"
            << "\t--------------- \n\n"
            << "\t\t t  : " << t_ << "\n"
            << "\t\t t' : " << t_prime_ << "\n"
            << "\t\t tz : " << t_perp_ << "\n"
            << "\t\t U  : " << U_ << "\n"
            << "\t\t U' : " << U_prime_ << "\n"
            << "\t\t V  : " << V_ << "\n"
            << "\t\t V' : " << V_prime_ << "\n"
            << "\n"
            << std::endl;
}
