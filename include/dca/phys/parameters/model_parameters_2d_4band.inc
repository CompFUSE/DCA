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
// Template specialization of model parameters for 2D 4-band model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::fourband_lattice<PointGroup>>> {
public:
  ModelParameters()
      : ei0_(0.),
        eb0_(0.),
        t0_(0.),
        ei1_(0.),
        eb1_(0.),
        t1_(0.),
        U0_(0.),
        U1_(0.),
        V_(0.),
        V_prime_(0.) {}

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
  double get_ei1() const {
    return ei1_;
  }
  double get_eb1() const {
    return eb1_;
  }
  double get_t1() const {
    return t1_;
  }
  double get_U0() const {
    return U0_;
  }
  double get_U1() const {
    return U1_;
  }
  double get_V() const {
    return V_;
  }
  double get_V_prime() const {
    return V_prime_;
  }

private:
  double ei0_;
  double eb0_;
  double t0_;
  double ei1_;
  double eb1_;
  double t1_;
  double U0_;
  double U1_;
  double V_;
  double V_prime_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::fourband_lattice<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(ei0_);
  buffer_size += concurrency.get_buffer_size(eb0_);
  buffer_size += concurrency.get_buffer_size(t0_);
  buffer_size += concurrency.get_buffer_size(ei1_);
  buffer_size += concurrency.get_buffer_size(eb1_);
  buffer_size += concurrency.get_buffer_size(t1_);
  buffer_size += concurrency.get_buffer_size(U0_);
  buffer_size += concurrency.get_buffer_size(U1_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::fourband_lattice<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, ei0_);
  concurrency.pack(buffer, buffer_size, position, eb0_);
  concurrency.pack(buffer, buffer_size, position, t0_);
  concurrency.pack(buffer, buffer_size, position, ei1_);
  concurrency.pack(buffer, buffer_size, position, eb1_);
  concurrency.pack(buffer, buffer_size, position, t1_);
  concurrency.pack(buffer, buffer_size, position, U0_);
  concurrency.pack(buffer, buffer_size, position, U1_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::fourband_lattice<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, ei0_);
  concurrency.unpack(buffer, buffer_size, position, eb0_);
  concurrency.unpack(buffer, buffer_size, position, t0_);
  concurrency.unpack(buffer, buffer_size, position, ei1_);
  concurrency.unpack(buffer, buffer_size, position, eb1_);
  concurrency.unpack(buffer, buffer_size, position, t1_);
  concurrency.unpack(buffer, buffer_size, position, U0_);
  concurrency.unpack(buffer, buffer_size, position, U1_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::fourband_lattice<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("fourband-model");

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
      reader_or_writer.execute("ei1", ei1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("eb1", eb1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t1", t1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U0", U0_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U1", U1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("V", V_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("V_prime", V_prime_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\nNo 4-band-model parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  std::cout << "\n\n"
            << "\t fourband-model : \n"
            << "\t--------------- \n\n"
            << "\t\t ei0  : " << ei0_ << "\n"
            << "\t\t eb0  : " << eb0_ << "\n"
            << "\t\t  t0  : " << t0_ << "\n"
            << "\t\t ei1  : " << ei1_ << "\n"
            << "\t\t eb1  : " << eb1_ << "\n"
            << "\t\t  t1  : " << t1_ << "\n"
            << "\t\t  U0  : " << U0_ << "\n"
            << "\t\t  U1  : " << U1_ << "\n"
            << "\t\t  V   : " << V_ << "\n"
            << "\t\t  V_prime : " << V_prime_ << "\n"
            << "\n"
            << std::endl;
}
