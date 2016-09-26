// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization of model parameters for tight-binding model.

template <typename Lattice>
class ModelParameters<models::TightBindingModel<Lattice>> {
public:
  ModelParameters() : t_(1.), t_prime_(0.), U_(4.), V_(0.), V_prime_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, int* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_t() const {
    return t_;
  }
  double get_t_prime() const {
    return t_prime_;
  }
  double get_U() const {
    return U_;
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
  double U_;
  double V_;
  double V_prime_;
};

template <typename Lattice>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<Lattice>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_);
  buffer_size += concurrency.get_buffer_size(t_prime_);
  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);

  return buffer_size;
}

template <typename Lattice>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<Lattice>>::pack(const Concurrency& concurrency,
                                                               int* buffer, int buffer_size,
                                                               int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_);
  concurrency.pack(buffer, buffer_size, position, t_prime_);
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
}

template <typename Lattice>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<Lattice>>::unpack(const Concurrency& concurrency,
                                                                 int* buffer, int buffer_size,
                                                                 int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_);
  concurrency.unpack(buffer, buffer_size, position, t_prime_);
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
}

template <typename Lattice>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<Lattice>>::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("2D-Hubbard-model");

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
    std::cout << "\nNo 2D-Hubbard-model parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}
