// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Template specialization of model parameters for the model described in "twoband_Cu.hpp".

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::TwoBandCu<PointGroup>>> {
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

  auto get_t00() const {
    return t00_;
  }
  auto get_tp00() const {
    return tp00_;
  }
  auto get_tpp00() const {
    return tpp00_;
  }

  auto get_t11() const {
    return t11_;
  }
  auto get_tp11() const {
    return tp11_;
  }
  auto get_tpp11() const {
    return tpp11_;
  }

  auto get_t01() const {
    return t01_;
  }
  auto get_tpp01() const {
    return tpp01_;
  }

  auto get_delta_E() const {
    return delta_E_;
  }

  auto get_U() const {
    return U_;
  }
  auto get_V() const {
    return V_;
  }
  auto get_J() const {
    return J_;
  }
  auto get_Jp() const {
    return Jp_;
  }

private:
  double t00_ = 0;
  double tp00_ = 0;
  double tpp00_ = 0;

  double t11_ = 0;
  double tp11_ = 0;
  double tpp11_ = 0;

  double t01_ = 0;
  double tpp01_ = 0;

  double delta_E_ = 0;

  double U_ = 0;
  double V_ = 0;

  double J_ = 0;
  double Jp_ = 0;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::TwoBandCu<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t00_);
  buffer_size += concurrency.get_buffer_size(tp00_);
  buffer_size += concurrency.get_buffer_size(tpp00_);
  buffer_size += concurrency.get_buffer_size(t11_);
  buffer_size += concurrency.get_buffer_size(tp11_);
  buffer_size += concurrency.get_buffer_size(tpp11_);
  buffer_size += concurrency.get_buffer_size(t01_);
  buffer_size += concurrency.get_buffer_size(tpp01_);
  buffer_size += concurrency.get_buffer_size(delta_E_);
  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(J_);
  buffer_size += concurrency.get_buffer_size(Jp_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::TwoBandCu<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t00_);
  concurrency.pack(buffer, buffer_size, position, tp00_);
  concurrency.pack(buffer, buffer_size, position, tpp00_);
  concurrency.pack(buffer, buffer_size, position, t11_);
  concurrency.pack(buffer, buffer_size, position, tp11_);
  concurrency.pack(buffer, buffer_size, position, tpp11_);
  concurrency.pack(buffer, buffer_size, position, t01_);
  concurrency.pack(buffer, buffer_size, position, tpp01_);
  concurrency.pack(buffer, buffer_size, position, delta_E_);
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, J_);
  concurrency.pack(buffer, buffer_size, position, Jp_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::TwoBandCu<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t00_);
  concurrency.unpack(buffer, buffer_size, position, tp00_);
  concurrency.unpack(buffer, buffer_size, position, tpp00_);
  concurrency.unpack(buffer, buffer_size, position, t11_);
  concurrency.unpack(buffer, buffer_size, position, tp11_);
  concurrency.unpack(buffer, buffer_size, position, tpp11_);
  concurrency.unpack(buffer, buffer_size, position, t01_);
  concurrency.unpack(buffer, buffer_size, position, tpp01_);
  concurrency.unpack(buffer, buffer_size, position, delta_E_);
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, J_);
  concurrency.unpack(buffer, buffer_size, position, Jp_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::TwoBandCu<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("TwoBandCu-model");

    auto try_to_read = [&](const std::string& name, auto& val) {
      try {
        reader_or_writer.execute(name, val);
      }
      catch (const std::exception& r_e) {
      }
    };

    try_to_read("t00", t00_);
    try_to_read("tp00", tp00_);
    try_to_read("tpp00", tpp00_);
    try_to_read("t11", t11_);
    try_to_read("tp11", tp11_);
    try_to_read("tpp11", tpp11_);
    try_to_read("t01", t01_);
    try_to_read("tpp01", tpp01_);

    try_to_read("delta-E", delta_E_);
    try_to_read("U", U_);
    try_to_read("V", V_);

    try_to_read("J", J_);
    try_to_read("Jp", Jp_);

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    throw(std::logic_error("No model parameters in input file."));
  }
}
