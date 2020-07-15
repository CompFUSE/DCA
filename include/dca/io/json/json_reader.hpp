// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Based on json-cpp.

#ifndef DCA_IO_JSON_JSON_READER_HPP
#define DCA_IO_JSON_JSON_READER_HPP

#include <cstdlib>
#include <complex>
#include <sstream>
#include <string>
#include <vector>
#include <type_traits>

#include "dca/util/type_utils.hpp"
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/json/json_parser/json_context.hpp"
#include "dca/io/json/json_parser/json_parser.hpp"
#include "dca/io/json/json_parser/whatever.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

namespace dca {
namespace io {
// dca::io::

class JSONReader {
public:
  typedef std::stringstream file_type;

  typedef Whatever JsonAccessor;
  typedef JSON_context JsonDataType;

public:
  JSONReader();

  bool is_reader() {
    return true;
  }
  bool is_writer() {
    return false;
  }

  std::string get_current_file_name() {
    return current_file_name;
  }

  void open_file(std::string file_name);
  void close_file() {}

  void open_group(std::string name) {
    my_paths.push_back(name);
  }
  void close_group() {
    my_paths.pop_back();
  }

  std::string get_path();

  template <typename arbitrary_struct_t>
  static void from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name);

  // Scalar type, std::string etc. 
  template <typename T>
  void execute(std::string name, T& value);

  // TODO: Remove? (only thing that depends on domains.hpp)
  template <typename domain_type>
  void execute(std::string /*name*/, func::dmn_0<domain_type>& /*dmn*/) {}

  template <typename scalartype, typename domain_type>
  void execute(func::function<scalartype, domain_type>& f);

  template <typename scalartype, typename domain_type>
  void execute(std::string name, func::function<scalartype, domain_type>& f);

  template <typename scalar_type>
  void execute(std::string name, dca::linalg::Vector<scalar_type, dca::linalg::CPU>& A);

  template <typename scalar_type>
  void execute(std::string name, dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A);

private:
  template <typename scalartype>
  void execute(std::string name, scalartype& value, const JsonAccessor& current_result,
               std::size_t index);

  template <typename scalartype, typename domain_type>
  void execute(std::string name, func::function<scalartype, domain_type>& f,
               const JsonAccessor& current_result, std::size_t index);

  template <typename scalartype, typename domain_type>
  void execute(std::string name, func::function<std::complex<scalartype>, domain_type>& f,
               const JsonAccessor& current_result, std::size_t index);

  template <typename scalar_type>
  void execute(std::string name, dca::linalg::Vector<scalar_type, dca::linalg::CPU>& A,
               const JsonAccessor& current_result, std::size_t index);

  template <typename scalar_type>
  void execute(std::string name, dca::linalg::Vector<std::complex<scalar_type>, dca::linalg::CPU>& A,
               const JsonAccessor& current_result, std::size_t index);

  template <typename scalar_type>
  void execute(std::string name, dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A,
               const JsonAccessor& current_result, std::size_t index);

  template <typename scalar_type>
  void execute(std::string name, dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU>& A,
               const JsonAccessor& current_result, std::size_t index);

private:
  void parse(std::string& file_name);

  const JsonAccessor& operator[](const std::string& key) const {
    return parse_result[key];
  }
  const JsonAccessor& operator[](int index) {
    return parse_result[index];
  }

  const JsonAccessor& get_JSON_tree() {
    return parse_result;
  }

private:
  std::string current_file_name;

  JSON_parser<JsonDataType> parser;
  const JsonAccessor& parse_result;

  std::vector<std::string> my_paths;
};

template <typename arbitrary_struct_t>
void JSONReader::from_file(arbitrary_struct_t& arbitrary_struct, std::string file_name) {
  JSONReader reader_obj;

  reader_obj.open_file(file_name);
  arbitrary_struct.read_write(reader_obj);
  reader_obj.close_file();
}

template <typename T>
void JSONReader::execute(std::string name, T& value) {
  std::size_t index = 0;
  if (index == my_paths.size())
    value <= parse_result[name];
  else {
    const JsonAccessor& new_result(parse_result[my_paths[index]]);

    execute(name, value, new_result, (++index));
  }
}

template <typename scalartype>
void JSONReader::execute(std::string name, scalartype& value, const JsonAccessor& current_result,
                         std::size_t index) {
  if (index == my_paths.size())
    value <= current_result[name];
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, value, new_result, (++index));
  }
}

template <typename scalartype, typename domain_type>
void JSONReader::execute(func::function<scalartype, domain_type>& f) {
  std::cout << "\t starts reading function : " << f.get_name() << "\n";
  execute(f.get_name(), f, parse_result, 0);
}

template <typename scalartype, typename domain_type>
void JSONReader::execute(std::string name, func::function<scalartype, domain_type>& f) {
  execute(name, f, parse_result, 0);
}

template <typename scalartype, typename domain_type>
void JSONReader::execute(std::string name, func::function<scalartype, domain_type>& f,
                         const JsonAccessor& current_result, std::size_t index) {
  if (index == my_paths.size()) {
    std::vector<std::vector<scalartype>> value;

    value <= current_result[name]["data"];

    int N_dmns = f.signature();
    for (std::size_t l = 0; l < value.size(); l++)
      f(l) = value[l][N_dmns];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, f, new_result, (++index));
  }
}

template <typename scalartype, typename domain_type>
void JSONReader::execute(std::string name, func::function<std::complex<scalartype>, domain_type>& f,
                         const JsonAccessor& current_result, std::size_t index) {
  const std::complex<scalartype> I(0, 1);

  if (index == my_paths.size()) {
    std::vector<std::vector<scalartype>> value;

    value <= current_result[name]["data"];

    int N_dmns = f.signature();
    for (std::size_t l = 0; l < value.size(); l++)
      f(l) = value[l][N_dmns] + I * value[l][N_dmns + 1];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, f, new_result, (++index));
  }
}

template <typename scalar_type>
void JSONReader::execute(std::string name, dca::linalg::Vector<scalar_type, dca::linalg::CPU>& V) {
  execute(name, V, parse_result, 0);
}

template <typename scalar_type>
void JSONReader::execute(std::string name, dca::linalg::Vector<scalar_type, dca::linalg::CPU>& V,
                         const JsonAccessor& current_result, std::size_t index) {
  if (index == my_paths.size()) {
    std::vector<scalar_type> value;

    value <= current_result[name]["data"];

    for (std::size_t i = 0; i < value.size(); i++)
      V[i] = value[i];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, V, new_result, (++index));
  }
}

template <typename scalar_type>
void JSONReader::execute(std::string name,
                         dca::linalg::Vector<std::complex<scalar_type>, dca::linalg::CPU>& V,
                         const JsonAccessor& current_result, std::size_t index) {
  if (index == my_paths.size()) {
    std::vector<scalar_type> value;

    value <= current_result[name]["data-real"];

    for (std::size_t i = 0; i < value.size(); i++)
      real(V[i]) = value[i];

    value <= current_result[name]["data-imag"];

    for (std::size_t i = 0; i < value.size(); i++)
      imag(V[i]) = value[i];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, V, new_result, (++index));
  }
}

template <typename scalar_type>
void JSONReader::execute(std::string name, dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A) {
  execute(name, A, parse_result, 0);
}

template <typename scalar_type>
void JSONReader::execute(std::string name, dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& A,
                         const JsonAccessor& current_result, std::size_t index) {
  if (index == my_paths.size()) {
    std::vector<std::vector<scalar_type>> value;

    value <= current_result[name]["data"];

    for (std::size_t i = 0; i < value.size(); i++)
      for (std::size_t j = 0; j < value[i].size(); j++)
        A(i, j) = value[i][j];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, A, new_result, (++index));
  }
}

template <typename scalar_type>
void JSONReader::execute(std::string name,
                         dca::linalg::Matrix<std::complex<scalar_type>, dca::linalg::CPU>& A,
                         const JsonAccessor& current_result, std::size_t index) {
  if (index == my_paths.size()) {
    std::vector<std::vector<scalar_type>> value;

    value <= current_result[name]["data-real"];

    for (std::size_t i = 0; i < value.size(); i++)
      for (std::size_t j = 0; j < value[i].size(); j++)
        real(A(i, j)) = value[i][j];

    value <= current_result[name]["data-imag"];

    for (std::size_t i = 0; i < value.size(); i++)
      for (std::size_t j = 0; j < value[i].size(); j++)
        imag(A(i, j)) = value[i][j];
  }
  else {
    const JsonAccessor& new_result(current_result[my_paths[index]]);

    execute(name, A, new_result, (++index));
  }
}

}  // io
}  // dca

#endif  // DCA_IO_JSON_JSON_READER_HPP
