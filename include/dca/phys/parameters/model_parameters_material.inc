// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization of model parameters for material models.

template <models::material_name_type name, typename PointGroup>
class ModelParameters<models::TightBindingModel<models::material_lattice<name, PointGroup>>> {
public:
  ModelParameters() : t_ij_file_name_("t_ij_file_name"), U_ij_file_name_("U_ij_file_name") {}

  template <typename Concurrency>
  int getBufferSize(Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(Concurrency& concurrency, int* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(Concurrency& concurrency, int* buffer, int buffer_size, int& position);

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
    Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_ij_file_name_);
  buffer_size += concurrency.get_buffer_size(U_ij_file_name_);

  return buffer_size;
}

template <models::material_name_type name, typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::material_lattice<name, PointGroup>>>::pack(
    Concurrency& concurrency, int* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_ij_file_name_);
  concurrency.pack(buffer, buffer_size, position, U_ij_file_name_);
}

template <models::material_name_type name, typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::material_lattice<name, PointGroup>>>::unpack(
    Concurrency& concurrency, int* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_ij_file_name_);
  concurrency.unpack(buffer, buffer_size, position, U_ij_file_name_);
}

template <models::material_name_type name, typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::material_lattice<name, PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("material-model");

    reader_or_writer.execute("t_ij-filename", t_ij_file_name_);
    reader_or_writer.execute("U_ij-filename", U_ij_file_name_);

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
    std::cout << "\nNo material-model parameters defined!\n" << std::endl;
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  std::cout << "\n\n"
            << "\tmaterial-model\n"
            << "\t--------------\n"
            << "\t\tt_ij-filename   : " << t_ij_file_name_ << "\n"
            << "\t\tU_ij-filename   : " << U_ij_file_name_ << "\n"
            << std::endl;
}
