// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Template specialization of model parameters for a 2D minimal model of iron oxypnictidesiron.
// See "Minimal two-band model of the superconducting iron oxypnictides,
// PHYSICAL REVIEW B 77, 220503 (2008)"

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::FeAsLattice<PointGroup>>> {
public:
  ModelParameters() = default;

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  // Tight binding parameters:
  // t1: nearest-neighbour, same orbital, in the orbital's \vec{l} direction.
  // t2: nearest-neighbour, same orbital, orthogonal direction respect to t1.
  // t3: next-nearest-neighbour, same orbital.
  // t4: next-nearest-neighbour, opposite orbital.
  double get_t1() const {
    return t1_;
  }
  double get_t2() const {
    return t2_;
  }
  double get_t3() const {
    return t3_;
  }
  double get_t4() const {
    return t4_;
  }

  // On-site, same orbital repulsion.
  double get_U() const {
    return U_;
  }
  // On-site, different orbital repulsion.
  double get_V() const {
    return V_;
  }
  // \vec{S_1}\vec{S_2} interaction.
  double get_J() const {
    return J_;
  }

  // Pair hoping interaction.
  double get_Jp() const {
    return Jp_;
  }

private:
  double t1_ = 0, t2_ = 0, t3_ = 0, t4_ = 0;
  double U_ = 0.;
  double V_ = 0.;
  double J_ = 0.;
  double Jp_ = 0.;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::FeAsLattice<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t1_);
  buffer_size += concurrency.get_buffer_size(t2_);
  buffer_size += concurrency.get_buffer_size(t3_);
  buffer_size += concurrency.get_buffer_size(t4_);
  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(J_);
  buffer_size += concurrency.get_buffer_size(Jp_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::FeAsLattice<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t1_);
  concurrency.pack(buffer, buffer_size, position, t2_);
  concurrency.pack(buffer, buffer_size, position, t3_);
  concurrency.pack(buffer, buffer_size, position, t4_);
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, J_);
  concurrency.pack(buffer, buffer_size, position, Jp_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::FeAsLattice<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t1_);
  concurrency.unpack(buffer, buffer_size, position, t2_);
  concurrency.unpack(buffer, buffer_size, position, t3_);
  concurrency.unpack(buffer, buffer_size, position, t4_);
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, J_);
  concurrency.unpack(buffer, buffer_size, position, Jp_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::FeAsLattice<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  auto tryToRead = [&](const std::string& name, auto& obj) {
    try {
      reader_or_writer.execute(name, obj);
    }
    catch (const std::exception& /*err*/) {
    }
  };

  try {
    reader_or_writer.open_group("FeAs-model");
    tryToRead("t1", t1_);
    tryToRead("t2", t2_);
    tryToRead("t3", t3_);
    tryToRead("t4", t4_);
    tryToRead("U", U_);
    tryToRead("V", V_);
    tryToRead("J", J_);
    tryToRead("Jp", Jp_);
    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\n  FeAs-model parameter group is not defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  if constexpr (ReaderOrWriter::is_reader)
    std::cout << "\n\n"
              << "\tFeAs-model : \n"
              << "\t--------------- \n\n"
              << "\t\t t  : " << t1_ << ", " << t2_ << ", " << t3_ << ", " << t4_ << "\n"
              << "\t\t U  : " << U_ << "\n"
              << "\t\t V  : " << V_ << "\n"
              << "\t\t J : " << J_ << "\n"
              << "\t\t Jp : " << Jp_ << "\n"
              << "\n"
              << std::endl;
}
