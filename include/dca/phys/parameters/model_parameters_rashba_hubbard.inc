// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Template specialization of model parameters for a 2D Rashba-Hubbard model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::RashbaHubbard<PointGroup>>> {
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

  double get_t() const {
    return t_;
  }

  double get_h() const {
    return h_;
  }

  double get_lambda() const {
    return lambda_;
  }

  double get_U() const {
    return U_;
  }

private:
  double t_ = 1;
  double h_ = 0.;
  double lambda_ = 0.;
  double U_ = 0.;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::RashbaHubbard<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_);
  buffer_size += concurrency.get_buffer_size(h_);
  buffer_size += concurrency.get_buffer_size(lambda_);
  buffer_size += concurrency.get_buffer_size(U_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::RashbaHubbard<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_);
  concurrency.pack(buffer, buffer_size, position, h_);
  concurrency.pack(buffer, buffer_size, position, lambda_);
  concurrency.pack(buffer, buffer_size, position, U_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::RashbaHubbard<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_);
  concurrency.unpack(buffer, buffer_size, position, h_);
  concurrency.unpack(buffer, buffer_size, position, lambda_);
  concurrency.unpack(buffer, buffer_size, position, U_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::RashbaHubbard<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  auto tryToRead = [&](const std::string& name, auto& obj) {
    try {
      reader_or_writer.execute(name, obj);
    }
    catch (const std::exception& /*err*/) {
      std::cerr << "Warning: Rashba-Hubbard-model has no " << name << " value specified.\n";
    }
  };

  try {
    reader_or_writer.open_group("Rashba-Hubbard-model");
    tryToRead("t", t_);
    tryToRead("h", h_);
    tryToRead("lambda", lambda_);
    tryToRead("U", U_);
    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\n  Rashba-Hubbard-model parameter group is not defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  if (ReaderOrWriter::is_reader) {
    std::cout << "\n\n"
              << "\tRasha-Hubbard-model : \n"
              << "\t--------------- \n\n"
              << "\t\t t  : " << t_ << "\n"
              << "\t\t h  : " << h_ << "\n"
              << "\t\t lambda  : " << lambda_ << "\n"
              << "\t\t U : " << U_ << "\n"
              << "\n"
              << std::endl;
  }
}


template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::RealRashbaHubbard<PointGroup>>> {
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

  double get_t() const {
    return t_;
  }

  double get_h() const {
    return h_;
  }

  double get_lambda() const {
    return lambda_;
  }

  double get_U() const {
    return U_;
  }

private:
  double t_ = 1;
  double h_ = 0.;
  double lambda_ = 0.;
  double U_ = 0.;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::RealRashbaHubbard<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_);
  buffer_size += concurrency.get_buffer_size(h_);
  buffer_size += concurrency.get_buffer_size(lambda_);
  buffer_size += concurrency.get_buffer_size(U_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::RealRashbaHubbard<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_);
  concurrency.pack(buffer, buffer_size, position, h_);
  concurrency.pack(buffer, buffer_size, position, lambda_);
  concurrency.pack(buffer, buffer_size, position, U_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::RealRashbaHubbard<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_);
  concurrency.unpack(buffer, buffer_size, position, h_);
  concurrency.unpack(buffer, buffer_size, position, lambda_);
  concurrency.unpack(buffer, buffer_size, position, U_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::RealRashbaHubbard<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  auto tryToRead = [&](const std::string& name, auto& obj) {
    try {
      reader_or_writer.execute(name, obj);
    }
    catch (const std::exception& /*err*/) {
      std::cerr << "Warning: Rashba-Hubbard-model has no " << name << " value specified.\n";
    }
  };

  try {
    reader_or_writer.open_group("Rashba-Hubbard-model");
    tryToRead("t", t_);
    tryToRead("h", h_);
    tryToRead("lambda", lambda_);
    tryToRead("U", U_);
    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\n  Rashba-Hubbard-model parameter group is not defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  if (ReaderOrWriter::is_reader) {
    std::cout << "\n\n"
              << "\tRasha-Hubbard-model : \n"
              << "\t--------------- \n\n"
              << "\t\t t  : " << t_ << "\n"
              << "\t\t h  : " << h_ << "\n"
              << "\t\t lambda  : " << lambda_ << "\n"
              << "\t\t U : " << U_ << "\n"
              << "\n"
              << std::endl;
  }
}
