// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization of model parameters for material models.

template <models::material_name_type name, typename PointGroup>
class ModelParameters<models::TightBindingModel<models::material_lattice<name, PointGroup>>> {
public:
  ModelParameters() : t_ij_file_name_("t_ij.txt"), U_ij_file_name_("U_ij.txt") {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  const std::string& get_t_ij_file_name() const {
    return t_ij_file_name_;
  }
  const std::string& get_U_ij_file_name() const {
    return U_ij_file_name_;
  }

  void set_t_ij_file_name(const std::string& file_name) {
    t_ij_file_name_ = file_name;
  }
  void set_U_ij_file_name(const std::string& file_name) {
    U_ij_file_name_ = file_name;
  }

private:
  std::string t_ij_file_name_;
  std::string U_ij_file_name_;
};

template <models::material_name_type name, typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::material_lattice<name, PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_ij_file_name_);
  buffer_size += concurrency.get_buffer_size(U_ij_file_name_);

  return buffer_size;
}

template <models::material_name_type name, typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::material_lattice<name, PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_ij_file_name_);
  concurrency.pack(buffer, buffer_size, position, U_ij_file_name_);
}

template <models::material_name_type name, typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::material_lattice<name, PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_ij_file_name_);
  concurrency.unpack(buffer, buffer_size, position, U_ij_file_name_);
}

template <models::material_name_type name, typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::material_lattice<name, PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("material-model");

    try {
      reader_or_writer.execute("t_ij-filename", t_ij_file_name_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U_ij-filename", U_ij_file_name_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}
